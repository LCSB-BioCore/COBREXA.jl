"""
    function max_min_driving_force(
        model::MetabolicModel,
        gibbs_free_energies::Dict{String,Float64},
        optimizer;
        ignore_metabolites::Vector{String} = [],
        constant_concentrations::Dict{String, Float64} = Dict{String, Float64}(),
        concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{Tuple{String,String},Float64}(),
        concentration_lb = 1e-6,
        concentration_ub = 10e-3,
        T = 298.15,
        R = 8.31446261815324e-3,
        modifications = [],
    )

Perform a max-min driving force analysis on the `model`, as defined by Noor,
Elad, et al. ("Pathway thermodynamics highlights kinetic obstacles in central
metabolism." *PLoS computational biology* 10.2, 2014).

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
thermodynamic calculations to assume constant pH and aqueous conditions.

`constant_concentrations` is used to fix the concentrations of certain
metabolites (such as the protons and CO₂). `concentration_ratios` is used to
specify additional constraints on metabolite pair concentrations (typically,
this is done with various cofactors such as for fixing the concentration ratios
of ADP and ATP). For example, you can optimistically fix the concentration of
ATP to be always 5× higher than of ADP by specifying
`Dict(("ATP","ADP") => 5.0)`

`concentration_lb` and `concentration_ub` set the `Cₗ` and `Cᵤ` in the
optimization problems.

`T` and `R` can be specified in the corresponding units; defaults are sensible
values in Kelvin and kJ/mol.
"""
function max_min_driving_force(
    model::MetabolicModel,
    gibbs_free_energies::Dict{String,Float64},
    optimizer;
    ignore_metabolites::Vector{String} = [],
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
    concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
        Tuple{String,String},
        Float64,
    }(),
    concentration_lb = 1e-6,
    concentration_ub = 10e-3,
    T = 298.15,
    R = 8.31446261815324e-3,
    modifications = [],
)
    rxn_index = Dict(rid => i for (i, rid) in enumerate(reactions(model)))

    # collect the information we need
    ridxs = [
        haskey(rxn_index, rid) ? rxn_index[rid] :
        throw(DomainError(rid, "reaction not found")) for
        rid in keys(gibbs_free_energies)
    ]
    dg0s = collect(values(gibbs_free_energies))

    # find the corresponding metabolites
    # (we use full S because we'll need it anyway)
    S = stoichiometry(model)
    midxs = Set{Int}()
    mets = metabolites(model)
    for ridx in ridxs
        for midx in findnz(S[:, ridx])[1]
            mets[midx] in ignore_metabolites || push!(midxs, midx)
        end
    end
    midxs = sort(collect(midxs))
    mids = mets[midxs]

    # reduce the stoichiometry to the reactions we're interested in
    S = S[midxs, ridxs]

    # start building the optimization model
    opt_model = Model(optimizer)
    @variables opt_model begin
        minDF
        logcs[1:length(midxs)]
        dgs[1:length(ridxs)]
    end

    RT = R * T
    @constraints opt_model begin
        minDF .<= -dgs
        dgs .<= 0
        dgs .== dg0s .+ RT .* S' * logcs
    end

    @objective(opt_model, Max, minDF)

    # add the absolute bounds
    all(in.(keys(constant_concentrations), Ref(mids))) || throw(
        DomainError(
            collect(keys(constant_concentrations)),
            "some of the metabolites not found in relevant reactions",
        ),
    )
    for (i, mid) in enumerate(mids)
        if haskey(constant_concentrations, mid)
            # we have an exact bound for this metabolite
            @constraint(opt_model, logcs[i] == log(constant_concentrations[mid]))
        else
            # this metabolite needs default bounds
            @constraint(
                opt_model,
                log(concentration_lb) <= logcs[i] <= log(concentration_ub)
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
        energies = Dict(
            rid => value(dgs[i]) for (i, rid) in enumerate(keys(gibbs_free_energies))
        ),
        concentrations = Dict(mid => exp(value(logcs[i])) for (i, mid) in enumerate(mids)),
    )
end
