"""
$(TYPEDSIGNATURES)

Build a max min driving force analysis model. See the docstring of
[`max_min_driving_force_analysis`](@ref) for details about the arguments.
"""
function build_max_min_driving_force_model(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    reference_flux = Dict{String,Float64}(),
    proton_ids::Vector{String},
    water_ids::Vector{String},
    concentration_lb = 1e-9, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
)
    # check if reactions that will be used in the model all have thermodynamic data, otherwise throw error.
    # only these reactions will form part of the model.
    rxn_source =
        isempty(reference_flux) ? A.reactions(model) : collect(keys(reference_flux))
    rxns = Set(x for x in rxn_source if !(x in ignore_reaction_ids))
    isempty(setdiff(rxns, keys(reaction_standard_gibbs_free_energies))) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            """
                      Not all reactions have thermodynamic data.
                      Either add the reactions with missing ΔG0s to ignore_reaction_ids,
                      or add the data to reaction_standard_gibbs_free_energies.
                      """,
        ),
    )
    rxns = collect(rxns)

    # import metabolite ids (if used), and reaction stoichiometries from an AbstractFBCModel
    mets = Set(m for rid in rxns for m in keys(A.reaction_stoichiometry(model, rid)))
    mets = collect(mets)
    stoi = A.stoichiometry(model)

    # create thermodynamic variables
    m = C.ConstraintTree(
        :max_min_driving_force^C.variable() +
        :log_metabolite_concentrations^C.variables(
            keys = Symbol.(mets),
            bounds = C.Between(log(concentration_lb), log(concentration_ub)),
        ) +
        :delta_G_reactions^C.variables(keys = Symbol.(rxns)),
    )

    #=
    Build gibbs free energy relations
    ΔG_rxn == ΔG0 + R * T * sum ( log_concentration_variable * stoichiometry_value )
    =#
    model_rxns = A.reactions(model)
    model_mets = A.metabolites(model)
    dG0s_met_ids_stoichs = Vector{Tuple{Float64,Vector{String},Vector{Float64}}}()
    for rxn in rxns # prepare to create CT below
        met_idxs, stoich_coeffs = findnz(stoi[:, findfirst(==(rxn), model_rxns)])
        met_ids = model_mets[met_idxs]
        dG0 = reaction_standard_gibbs_free_energies[rxn]
        push!(dG0s_met_ids_stoichs, (dG0, met_ids, stoich_coeffs))
    end

    m *=
        :delta_G_reaction_equations^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = -m.delta_G_reactions[Symbol(rxn)].value +
                        dG0 +
                        R *
                        T *
                        sum(
                            m.log_metabolite_concentrations[Symbol(met_id)].value * stoich for (met_id, stoich) in zip(met_ids, stoich_coeffs)
                        ),
                bound = C.EqualTo(0.0),
            ) for (rxn, (dG0, met_ids, stoich_coeffs)) in zip(rxns, dG0s_met_ids_stoichs)
        )

    #=
    Set proton log concentration to zero so that it won't impact any
    calculations (biothermodynamics assumption). Also set water concentrations
    to zero (aqueous condition assumptions). How true is "aqueous conditions"?
    Debatable...
    =#
    for met in [Symbol.(proton_ids); Symbol.(water_ids)]
        if haskey(m.log_metabolite_concentrations, met)
            m.log_metabolite_concentrations[met] = C.Constraint(
                value = m.log_metabolite_concentrations[met].value,
                bound = C.EqualTo(0.0),
            )
        end
    end

    #=
    Add thermodynamic feasibility constraint (ΔG < 0 for a feasible reaction in flux direction).
    Add objective constraint to solve max min problem.
    =#
    m *=
        :reaction_delta_G_margin^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = m.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(reference_flux, rxn, 1.0)),
                bound = C.Between(-Inf, 0.0),
            ) for rxn in rxns
        )

    m *=
        :min_driving_force_margin^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = m.max_min_driving_force.value +
                        m.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(reference_flux, rxn, 1.0)),
                bound = C.Between(-Inf, 0.0),
            ) for rxn in rxns
        )

    m
end

export build_max_min_driving_force_model


"""
$(TYPEDSIGNATURES)

Add constraints to `m` that represents ratios of variables in log space:
`log(x/y) = log(const)` where `x` and `y` are variables specified by `on`. These
constraints are called `name` in the resultant model. The constraints are
specified by `ratios`, which is a dictionary mapping a constraint id to a tuple
which consists of the variable ids, `x`, `y`, and the ratio value, `const`. The
latter is logged internally.

# Example
```
m *= :log_ratio_constraints^log_ratio_constraints(
    Dict("atp" => ("atp_c", "adp_c", log(10.0)), "nadh" => ("nadh_c", "nad_c", log(0.13))),
    m.log_metabolite_concentrations,
)
```
"""
log_ratio_constraints(
    ratios::Dict{String,Tuple{String,String,Float64}},
    on::C.ConstraintTree,
) = C.ConstraintTree(
    Symbol(cid) => C.Constraint(
        value = on[Symbol(var1)].value - on[Symbol(var2)].value,
        bound = C.EqualTo(log(ratio)),
    ) for (cid, (var1, var2, ratio)) in ratios
)

export log_ratio_constraints

"""
$(TYPEDSIGNATURES)

Perform a max-min driving force analysis using `optimizer` on the `model` with
supplied ΔG⁰s in `reaction_standard_gibbs_free_energies`, as defined by Noor, et
al., "Pathway thermodynamics highlights kinetic obstacles in central
metabolism.", PLoS computational biology, 2014.

Optionally, `reference_flux` can be used to set the directions of each reaction
in `model` (all reactions are assumed to proceed forward by default). The
supplied `reference_flux` should be free of internal cycles i.e.
thermodynamically consistent. This optional input is important if a reaction in
`model` normally runs in reverse (negative flux). Note, only the signs are
extracted from this input, so be careful with floating point precision near 0.

The max-min driving force algorithm returns the Gibbs free energy of the
reactions, the concentrations of metabolites and the actual maximum minimum
driving force. The optimization problem solved is:
```
max min -ΔᵣG
s.t. ΔᵣG = ΔrG⁰ + R T S' ln(C)
     ΔᵣG ≤ 0
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
where `ΔrG` are the Gibbs energies dissipated by the reactions, R is the gas
constant, T is the temperature, S is the stoichiometry of the model, and C is
the vector of metabolite concentrations (and their respective lower and upper
bounds).

In case no feasible solution exists, `nothing` is returned.

Reactions specified in `ignore_reaction_ids` are internally ignored when
calculating the max-min driving force. Importantly, this should include water
and proton importers.

Since biochemical thermodynamics are assumed, the `proton_ids` and `water_ids`
need to be specified so that they can be ignored in the calculations.
Effectively this assumes an aqueous environment at constant pH is used.

`constant_concentrations` is used to fix the concentrations of certain
metabolites (such as CO₂) by passing a dictionary mapping metabolite id to its
constant value. `concentration_ratios` is used to specify additional constraints
on metabolite pair concentrations (typically, this is done with various
cofactors, such as the ATP/ADP ratio. For example, you can fix the concentration
of ATP to be always 5× higher than of ADP by specifying `Dict("atp_ratio" =>
("atp_c","adp_c", 5.0))` if these metabolites are called `atp_c` and `adp_c` in
the model. `concentration_lb` and `concentration_ub` set the `Cₗ` and `Cᵤ` in
the optimization problems (these are overwritten by `constant_concentrations` if
supplied).

`T` and `R` can be specified in the corresponding units; defaults are K and
kJ/K/mol. The unit of metabolite concentrations is typically molar. As usual,
optimizer settings can be changed with `modifications`.
"""
function max_min_driving_force_analysis(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    reference_flux = Dict{String,Float64}(),
    concentration_ratios = Dict{String,Tuple{String,String,Float64}}(),
    constant_concentrations = Dict{String,Float64}(),
    proton_ids,
    water_ids,
    concentration_lb = 1e-9, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
    modifications = [],
    optimizer,
)
    m = build_max_min_driving_force_model(
        model,
        reaction_standard_gibbs_free_energies;
        reference_flux,
        concentration_lb,
        concentration_ub,
        R,
        T,
        ignore_reaction_ids,
        water_ids,
        proton_ids,
    )

    for (mid, val) in constant_concentrations
        m.log_metabolite_concentrations[Symbol(mid)].bound = C.EqualTo(log(val))
    end

    m *=
        :metabolite_ratio_constraints^log_ratio_constraints(
            concentration_ratios,
            m.log_metabolite_concentrations,
        )


    optimized_constraints(
        m;
        objective = m.max_min_driving_force.value,
        optimizer,
        modifications,
    )
end

export max_min_driving_force_analysis
