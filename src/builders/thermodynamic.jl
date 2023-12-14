"""
The function uses the supplied `optimizer` and `reaction_standard_gibbs_free_energies`.
Optionally, `flux_solution` can be used to set the directions of each reaction in `model`
(all reactions are assumed to proceed forward and are active by default). The supplied
`flux_solution` should be free of internal cycles i.e. thermodynamically consistent. This
optional input is important if a reaction in `model` normally runs in reverse (negative
flux). Note, reactions in `flux_solution` that are smaller than `small_flux_tol` are also
ignored in the analysis function (for numerical stability).


Reactions specified in `ignore_reaction_ids` are internally ignored when calculating the
max-min driving force. This should include water and proton importers.

Since biochemical thermodynamics are assumed, the `proton_ids` and `water_ids` need to be
specified so that they can be ignored in the calculations. Effectively this assumes an
aqueous environment at constant pH is used.

`constant_concentrations` is used to fix the concentrations of certain metabolites (such as
CO₂). `concentration_ratios` is used to specify additional constraints on metabolite pair
concentrations (typically, this is done with various cofactors such as the ATP/ADP ratio.
For example, you can fix the concentration of ATP to be always 5× higher than of ADP by
specifying `Dict(("ATP","ADP") => 5.0)`

`concentration_lb` and `concentration_ub` set the `Cₗ` and `Cᵤ` in the
optimization problems.

`T` and `R` can be specified in the corresponding units; defaults are K and kJ/K/mol.
"""    
function max_min_driving_force(
    model::MetabolicModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64},
    optimizer;
    flux_solution::Dict{String,Float64} = Dict{String,Float64}(),
    proton_ids::Vector{String} = ["h_c", "h_e"],
    water_ids::Vector{String} = ["h2o_c", "h2o_e"],
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
    concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
        Tuple{String,String},
        Float64,
    }(),
    concentration_lb = 1e-9,
    concentration_ub = 100e-3,
    T = _constants.T,
    R = _constants.R,
    small_flux_tol = 1e-6,
    modifications = [],
    ignore_reaction_ids = [],
)
    opt_model = Model(optimizer)

    @variables opt_model begin
        mmdf
        logcs[1:n_metabolites(model)]
        dgrs[1:n_reactions(model)]
    end

    # set proton log concentration to zero so that it won't impact any calculations (biothermodynamics assumption)
    proton_idxs = Int.(indexin(proton_ids, metabolites(model)))
    for idx in proton_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # set water concentrations to zero (aqueous condition assumptions)
    water_idxs = Int.(indexin(water_ids, metabolites(model)))
    for idx in water_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # only consider reactions with supplied thermodynamic data AND have a flux bigger than
    # small_flux_tol => finds a thermodynamic profile that explains flux_solution
    active_rids = filter(
        rid ->
            haskey(reaction_standard_gibbs_free_energies, rid) &&
                abs(get(flux_solution, rid, small_flux_tol / 2)) > small_flux_tol &&
                !(rid in ignore_reaction_ids),
        reactions(model),
    )
    active_ridxs = Int.(indexin(active_rids, reactions(model)))

    # give dummy dG0 for reactions that don't have data
    dg0s =
        [get(reaction_standard_gibbs_free_energies, rid, 0.0) for rid in reactions(model)]

    S = stoichiometry(model)

    @constraint(opt_model, dgrs .== dg0s .+ (R * T) * S' * logcs)

    # thermodynamics should correspond to the fluxes
    flux_signs = [sign(get(flux_solution, rid, 1.0)) for rid in reactions(model)]

    # only constrain reactions that have thermo data
    @constraint(opt_model, dgrs[active_ridxs] .* flux_signs[active_ridxs] .<= 0)

    # add the absolute bounds
    missing_mets =
        [mid for mid in keys(constant_concentrations) if !(mid in metabolites(model))]
    !isempty(missing_mets) &&
        throw(DomainError(missing_mets, "metabolite(s) not found in model."))
    for (midx, mid) in enumerate(metabolites(model)) # idx in opt_model (missing ignore_metabolites)
        midx in water_idxs && continue
        midx in proton_idxs && continue
        if haskey(constant_concentrations, mid)
            JuMP.fix(logcs[midx], log(constant_concentrations[mid]))
        else
            # this metabolite needs bounds
            @constraint(
                opt_model,
                log(concentration_lb) <= logcs[midx] <= log(concentration_ub)
            )
        end
    end

    # add the relative bounds
    for ((mid1, mid2), val) in concentration_ratios
        idxs = indexin([mid1, mid2], metabolites(model)) # TODO: this is not performant
        any(isnothing.(idxs)) &&
            throw(DomainError((mid1, mid2), "metabolite pair not found in model."))
        @constraint(opt_model, logcs[idxs[1]] == log(val) + logcs[idxs[2]])
    end

    @constraint(opt_model, mmdf .<= -dgrs[active_ridxs] .* flux_signs[active_ridxs])

    @objective(opt_model, Max, mmdf)

    # apply the modifications, if any
    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)

    is_solved(opt_model) || return nothing

    return (
        mmdf = value(opt_model[:mmdf]),
        dg_reactions = Dict(
            rid => value(opt_model[:dgrs][i]) for (i, rid) in enumerate(reactions(model))
        ),
        concentrations = Dict(
            mid => exp(value(opt_model[:logcs][i])) for
            (i, mid) in enumerate(metabolites(model))
        ),
    )
end