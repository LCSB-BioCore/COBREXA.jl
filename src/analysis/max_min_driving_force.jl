"""
    max_min_driving_force(
        model::MetabolicModel,
        gibbs_free_energies::Dict{String,Float64},
        optimizer;
        ignore_metabolites::Vector{String} = ["h", "h2o"],
        constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
        concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
            Tuple{String,String},
            Float64,
        }(),
        concentration_lb = 1e-6,
        concentration_ub = 10e-3,
        T = _constants.T,
        R = _constants.R,
        modifications = [],
    )

Perform a max-min driving force analysis on the `model`, as defined by Noor, et al.,
"Pathway thermodynamics highlights kinetic obstacles in central metabolism.", PLoS
computational biology, 2014.

The analysis uses the supplied `optimizer` and Gibbs free energies of the
reactions (in `gibbs_free_energies`) to find the max-min driving force, Gibbs
free energy of the reactions and the concentrations of metabolites that
optimize the following problem:
```
max min -ΔᵣG
s.t. ΔᵣG = ΔᵣG⁰ + R T S' ln(C)
     ΔᵣG ≤ 0 (∀r)
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
where `ΔᵣG` are the Gibbs energies dissipated by the reactions, `ΔᵣG⁰` are the
Gibbs free energies of the reactions, R is the gas constant, T is the
temperature, S is the stoichiometry of the model, and C is the vector of
metabolite concentrations (and their respective lower and upper bounds).

In case no feasible solution exists, `nothing` is returned.

Metabolites specified in `ignore_metabolites` are internally ignored -- that
allows to specify e.g. removal of protons and water, thus allowing the
thermodynamic calculations to assume constant pH and aqueous conditions. Note, if using
biochemical thermodynamic data then you _must_ include the ids of protons and water here.

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
    gibbs_free_energies::Dict{String,Float64},
    optimizer;
    ignore_metabolites::Vector{String} = ["h", "h2o"],
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
    concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
        Tuple{String,String},
        Float64,
    }(),
    concentration_lb = 1e-6,
    concentration_ub = 10e-3,
    T = _constants.T,
    R = _constants.R,
    modifications = [],
)
    # get ΔGs in order of reactions
    dg0s = [gibbs_free_energies[rid] for rid in reactions(model)]

    # get all metabolites that participate in the ΔG calculations (no protons, no water)
    midxs = [
        midx for
        (midx, mid) in enumerate(metabolites(model)) if !(mid in ignore_metabolites)
    ] # idx in original model
    mids = metabolites(model)[midxs]

    # find the corresponding metabolites
    S = stoichiometry(model)[midxs, :]

    # start building the optimization model
    opt_model = Model(optimizer)
    @variables opt_model begin
        minDF
        logcs[1:length(midxs)]
        dgs[1:n_reactions(model)]
    end

    @constraints opt_model begin
        minDF .<= -dgs
        dgs .<= 0
        dgs .== dg0s .+ R * T .* S' * logcs
    end

    @objective(opt_model, Max, minDF)

    # add the absolute bounds
    for (midx, mid) in enumerate(mids) # idx in opt_model (missing ignore_metabolites)
        if haskey(constant_concentrations, mid)
            # we have an exact bound for this metabolite
            @constraint(opt_model, logcs[midx] == log(constant_concentrations[mid]))
        else
            # this metabolite needs default bounds
            @constraint(
                opt_model,
                log(concentration_lb) <= logcs[midx] <= log(concentration_ub)
            )
        end
    end

    # add the relative bounds
    for ((mid1, mid2), val) in concentration_ratios
        idxs = indexin([mid1, mid2], mids)
        any(isnothing.(idxs)) && throw(
            DomainError((mid1, mid2), "metabolite pair not found in relevant reactions"),
        )
        @constraint(opt_model, logcs[idxs[1]] == log(val) + logcs[idxs[2]])
    end

    # apply the modifications, if any
    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)

    is_solved(opt_model) || return nothing

    return (
        mmdf = objective_value(opt_model),
        dgs = Dict(rid => value(dgs[i]) for (i, rid) in enumerate(reactions(model))),
        cs = Dict(mid => exp(value(logcs[i])) for (i, mid) in enumerate(mids)),
    )
end
