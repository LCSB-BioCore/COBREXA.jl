"""
    max_min_driving_force(
        model::StandardModel,
        optimizer,
        thermodynamic_data;
        modifications=[],
        proton_id="h_c",
        water_id="h2o_c",
        concentration_ratios=[("atp_c", "adp_c", 10.0),
        ("adp_c", "amp_c", 1.0),
        ("nadph_c", "nadp_c", 10.0),
        ("nadh_c", "nad_c", 0.1)],
        constant_concentrations=[("coa_c", 1e-3),
        ("co2_c", 10e-6),
        ("pi_c", 10e-3),
        ("ppi_c", 1e-3)],
        concentration_lb=1e-6,
        concentration_ub=10e-3)

Perform max min driving force analysis on `model` using `optimizer` and
`thermodynamic_data`. The `optimizer` can be modified using `modifications`. Returns the
maximum minimum driving force, the Gibbs free energy of reactions and the concentrations
of metabolites that solve
```
max min -ΔᵣG'
s.t. ΔᵣG' = ΔᵣG'⁰ + R*T*S'*ln(C)
     ΔᵣG' <= 0 ∀ r
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
See `Noor, Elad, et al. "Pathway thermodynamics highlights kinetic obstacles in central
metabolism." PLoS computational biology 10.2 (2014): e1003483.` for more information.

Note, protons and water need to be removed from the analysis because they do not figure into
the thermodynamic calculations (constant pH and aqueous conditions are assumed). Typically,
cofactors such as ATP, ADP, etc. are constrained by their ratios, as in
`concentration_ratios`, which is a vector of tuples like `[(numerator, denominator,
value),...]`. For the first element this corresponds to `numerator/denominator = value`.
Alternatively, metabolites in `constant_concentrations` can be directly constrained to
specific values, with the format being a vector of tuples `[(metabolite, concentration),...]`.
Sensible defaults are supplied here, although the name space needs to be updated depending
on the model. Finally, `Cₗ` and `Cᵤ` are set with `concentration_lb` and `concentration_ub`.

# Example
```
mmdf, dgs, concens = max_min_driving_force(
    model,
    Tulip.Optimizer,
    thermodynamic_data;
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
    optimizer,
    thermodynamic_data;
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
    rids = filter(x -> haskey(thermodynamic_data, x), reactions(model))
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

    dg0s = [thermodynamic_data[rid] for rid in rids]

    # modify optimization problem
    opt_model = Model(optimizer)
    for mod in modifications
        mod(model, opt_model)
    end

    @variables opt_model begin
        minDF
        dgs[1:length(dg0s)]
        logcs[1:size(S, 2)]
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
            @constraint(opt_model, logcs[first(i)] == log(constant_concentration_dict[mid]))
        else
            i == indexin([mid], mids)
            isnothing(i) && continue
            @constraint(opt_model, log_lb <= logcs[i] <= log_ub)
        end
    end

    for (met1, met2, val) in concentration_ratios
        i, j = indexin([met1, met2], mids)
        isnothing(i) || isnothing(j) && continue
        @constraint(opt_model, logcs[first(i)] == log(val) + logcs[first(j)])
    end

    @objective(opt_model, Max, minDF)

    optimize!(opt_model)

    is_solved(opt_model) || return nothing, nothing, nothing

    return objective_value(opt_model),
    Dict(rid => value(dgs[i]) for (i, rid) in enumerate(rids)),
    Dict(mid => exp(value(logcs[i])) for (i, mid) in enumerate(mids))
end
