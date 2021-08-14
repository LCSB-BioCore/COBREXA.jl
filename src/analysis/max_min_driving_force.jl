"""
    max_min_driving_force(
        model::StandardModel,
        standard_gibbs_reaction_energies::Dict{String, Float64},
        optimizer;
        modifications = [],
        proton_id = "h_c",
        water_id = "h2o_c",
        concentration_ratios = [
            ("atp_c", "adp_c", 10.0),
            ("adp_c", "amp_c", 1.0),
            ("nadph_c", "nadp_c", 10.0),
            ("nadh_c", "nad_c", 0.1),
        ],
        constant_concentrations = [
            ("coa_c", 1e-3),
            ("co2_c", 10e-6),
            ("pi_c", 10e-3),
            ("ppi_c", 1e-3),
        ],
        concentration_lb = 1e-6,
        concentration_ub = 10e-3,
)

Perform max min driving force analysis on `model` using `optimizer` and
`standard_gibbs_reaction_energies`. The `optimizer` can be modified using `modifications` but not
underlying optimization problem. Returns the maximum minimum driving force, the Gibbs free
energy of reactions and the concentrations of metabolites that solve
```
max min -ΔᵣG'
s.t. ΔᵣG' = ΔᵣG'⁰ + R*T*S'*ln(C)
     ΔᵣG' <= 0 ∀ r
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
See `Noor, Elad, et al. "Pathway thermodynamics highlights kinetic obstacles in central
metabolism." PLoS computational biology 10.2 (2014): e1003483.` for more information.

Internally protons and water are removed from the model because biological thermodynamic
calculations assume constant pH and aqueous conditions. Typically, cofactors such as ATP,
ADP, etc. are constrained by their ratios, as in `concentration_ratios`, which is a vector
of tuples like `[(numerator, denominator, value),...]`. For the first element this
corresponds to `numerator/denominator = value`. Alternatively, metabolites in
`constant_concentrations` can be directly constrained to specific values, with the format
being a vector of tuples `[(metabolite, concentration),...]`. Sensible defaults are supplied
here, although the name space needs to be updated depending on the model. Finally, `Cₗ` and
`Cᵤ` are set with `concentration_lb` and `concentration_ub`.

# Example
```
mmdf, dgs, concens = max_min_driving_force(
    model,
    gibbs_free_energies,
    Tulip.Optimizer;
    proton_id = "h",
    water_id = "h2o",
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 500)],
    concentration_ratios = [("atp", "adp", 10.0), ("nadh", "nad", 0.1)],
    constant_concentrations = [("pi", 10e-3)],
    concentration_lb = 1e-6,
    concentration_ub = 10e-3,
)
```
"""
function max_min_driving_force(
    model::StandardModel,
    standard_gibbs_reaction_energies::Dict{String,Float64},
    optimizer;
    modifications = [],
    proton_id = "h_c",
    water_id = "h2o_c",
    concentration_ratios = [
        ("atp_c", "adp_c", 10.0),
        ("adp_c", "amp_c", 1.0),
        ("nadph_c", "nadp_c", 10.0),
        ("nadh_c", "nad_c", 0.1),
    ],
    constant_concentrations = [
        ("coa_c", 1e-3),
        ("co2_c", 10e-6),
        ("pi_c", 10e-3),
        ("ppi_c", 1e-3),
    ],
    concentration_lb = 1e-6,
    concentration_ub = 10e-3,
)

    # find reactions with thermodynamic data, ignore all other reactions in model
    rids = filter(x -> haskey(standard_gibbs_reaction_energies, x), reactions(model))
    ridxs = Int.(indexin(rids, reactions(model)))

    # remove protons, water and all metabolites not involved in reactions that have thermodynamic data
    mids = unique(
        vcat(
            [
                collect(keys(rxn.metabolites)) for
                (rid, rxn) in model.reactions if rid in rids
            ]...,
        ),
    )
    filter!(x -> !(x in [proton_id, water_id]), mids)
    midxs = Int.(indexin(mids, metabolites(model)))

    St = (stoichiometry(model)[midxs, ridxs])'

    RT = 298.15 * 8.314e-3 # kJ/mol

    dg0s = [standard_gibbs_reaction_energies[rid] for rid in rids]

    # modify optimization problem
    opt_model = Model(optimizer)
    for mod in modifications
        mod(model, opt_model)
    end

    @variables opt_model begin
        minDF
        dgs[1:length(dg0s)]
        logcs[1:size(St, 2)]
    end

    @constraints opt_model begin
        minDF .<= -dgs
        dgs .<= 0
        dgs .== dg0s .+ RT .* St * logcs
    end

    log_lb = log(concentration_lb)
    log_ub = log(concentration_ub)
    constant_concentration_dict = Dict(met => val for (met, val) in constant_concentrations)
    for (i, mid) in enumerate(mids)
        if haskey(constant_concentration_dict, mid)
            @constraint(opt_model, logcs[i] == log(constant_concentration_dict[mid]))
        else
            @constraint(opt_model, log_lb <= logcs[i] <= log_ub)
        end
    end

    for (met1, met2, val) in concentration_ratios
        i = first(indexin([met1], mids))
        isnothing(i) && continue
        j = first(indexin([met2], mids))
        isnothing(j) && continue
        @constraint(opt_model, logcs[i] == log(val) + logcs[j])
    end

    @objective(opt_model, Max, minDF)

    optimize!(opt_model)

    is_solved(opt_model) || return nothing

    return (
        max_min_driving_force = objective_value(opt_model),
        optimal_gibbs_free_energies = Dict(
            rid => value(dgs[i]) for (i, rid) in enumerate(rids)
        ),
        optimal_concentrations = Dict(
            mid => exp(value(logcs[i])) for (i, mid) in enumerate(mids)
        ),
    )
end
