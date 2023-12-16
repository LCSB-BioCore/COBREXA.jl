"""
$(TYPEDSIGNATURES)

The function uses the supplied `optimizer` and
`reaction_standard_gibbs_free_energies`. Optionally, `flux_solution` can be used
to set the directions of each reaction in `model` (all reactions are assumed to
proceed forward and are active by default). The supplied `flux_solution` should
be free of internal cycles i.e. thermodynamically consistent. This optional
input is important if a reaction in `model` normally runs in reverse (negative
flux). Note, reactions in `flux_solution` that are smaller than `small_flux_tol`
are also ignored in the analysis function (for numerical stability).


Reactions specified in `ignore_reaction_ids` are internally ignored when
calculating the max-min driving force. This should include water and proton
importers.

Since biochemical thermodynamics are assumed, the `proton_ids` and `water_ids`
need to be specified so that they can be ignored in the calculations.
Effectively this assumes an aqueous environment at constant pH is used.

`constant_concentrations` is used to fix the concentrations of certain
metabolites (such as CO₂). `concentration_ratios` is used to specify additional
constraints on metabolite pair concentrations (typically, this is done with
various cofactors such as the ATP/ADP ratio. For example, you can fix the
concentration of ATP to be always 5× higher than of ADP by specifying
`Dict(("ATP","ADP") => 5.0)`

`concentration_lb` and `concentration_ub` set default concentration bounds, in M
by default.

`T` and `R` can be specified in the corresponding units; defaults are K and
kJ/K/mol.
"""
function build_max_min_driving_force_model(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    flux_solution = Dict{String,Float64}(),
    proton_ids = ["h_c", "h_e"],
    water_ids = ["h2o_c", "h2o_e"],
    concentration_lb = 1e-9,
    concentration_ub = 100e-3,
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
)
    # check if reactions that will be used in the model all have thermodynamic data, otherwise throw error.
    # only these reactions will form part of the model.
    rxns = filter(
        x -> !(x in ignore_reaction_ids),
        isempty(flux_solution) ? A.reactions(model) : collect(keys(flux_solution)),
    )
    all(in.(rxns, Ref(collect(keys(reaction_standard_gibbs_free_energies))))) || throw(
        ArgumentError("""
                      Not all reactions have thermodynamic data. 
                      Either add the reactions with missing ΔG0s to ignore_reaction_ids, 
                      or add the data to reaction_standard_gibbs_free_energies.
                      """),
    )

    # import metabolite ids (if used), and reaction stoichiometries from an AbstractFBCModel
    mets = Set(Symbol[])
    for rid in rxns
        for met in keys(A.reaction_stoichiometry(model, rid))
            push!(mets, Symbol(met))
        end
    end
    mets = [m for m in mets] # TODO: constrainttrees #12
    stoi = A.stoichiometry(model)

    # create thermodynamic variables
    m = C.ConstraintTree(:max_min_driving_force^C.variable())
    m +=
        :log_metabolite_concentrations^C.variables(
            keys = mets,
            bounds = zip(log(concentration_lb), log(concentration_ub)),
        )
    m += :delta_G_reactions^C.variables(keys = Symbol.(rxns))

    #= 
    Build gibbs free energy relations
    ΔG_rxn == ΔG0 + R * T * sum ( log_concentration_variable * stoichiometry_value )
    =#
    model_rxns = A.reactions(model)
    model_mets = A.metabolites(model)
    for rxn in rxns
        met_idxs, stoich_coeffs = findnz(stoi[:, findfirst(==(rxn), model_rxns)])
        met_ids = model_mets[met_idxs]

        dG0 = reaction_standard_gibbs_free_energies[rxn]

        m *=
            :delta_G_reaction_equations^Symbol(
                rxn,
            )^C.Constraint(
                value = -m.delta_G_reactions[Symbol(rxn)].value +
                        dG0 +
                        R *
                        T *
                        sum(
                            m.log_metabolite_concentrations[Symbol(met_id)].value *
                            stoich for (met_id, stoich) in zip(met_ids, stoich_coeffs)
                        ),
                bound = 0.0,
            )
    end

    #= 
    Set proton log concentration to zero so that it won't impact any
    calculations (biothermodynamics assumption). Also set water concentrations
    to zero (aqueous condition assumptions). How true is "aqueous conditions"? 
    Debatable... 
    =#
    log_met_conc_zero = [Symbol.(proton_ids); Symbol.(water_ids)]
    for met in keys(m.log_metabolite_concentrations)
        if met in log_met_conc_zero
            m.log_metabolite_concentrations[met].bound = 0.0
        end
    end

    #= 
    Add thermodynamic feasibility constraint (ΔG < 0 for a feasible reaction in flux direction).
    Add objective constraint to solve max min problem.
    =#
    for rxn in rxns
        m *=
            :thermodynamic_feasibility^Symbol(
                rxn,
            )^C.Constraint(
                value = m.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(flux_solution, rxn, 1.0)),
                bound = (-Inf, 0.0),
            )

        m *=
            :max_min_constraints^Symbol(
                rxn,
            )^C.Constraint(
                value = m.max_min_driving_force.value +
                        m.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(flux_solution, rxn, 1.0)),
                bound = (-Inf, 0.0),
            )
    end

    m
end

export build_max_min_driving_force_model

function add_metabolite_ratio_constraints!(
    m::C.ConstraintTree,
    concentration_ratios::Dict{String,Tuple{String,String,Float64}},
)
    for (constraint_id, (met1, met2, ratio)) in concentration_ratios
        m *=
            :metabolite_ratio_constraints^Symbol(
                constraint_id,
            )^C.Constraint(
                value = m.log_metabolite_concentrations[Symbol(met1)].value -
                        m.log_metabolite_concentrations[Symbol(met2)].value,
                bound = log(ratio),
            )
    end

    m
end

export add_metabolite_ratio_constraints!

"""
$(TYPEDSIGNATURES)

Perform a max-min driving force analysis on the `model`, as defined by Noor, et al.,
"Pathway thermodynamics highlights kinetic obstacles in central metabolism.", PLoS
computational biology, 2014.

The function uses the supplied `optimizer` and `reaction_standard_gibbs_free_energies`.
Optionally, `flux_solution` can be used to set the directions of each reaction in `model`
(all reactions are assumed to proceed forward and are active by default). The supplied
`flux_solution` should be free of internal cycles i.e. thermodynamically consistent. This
optional input is important if a reaction in `model` normally runs in reverse (negative
flux). Note, reactions in `flux_solution` that are smaller than `small_flux_tol` are also
ignored in the analysis function (for numerical stability).

The max-min driving force algorithm returns the Gibbs free energy of the reactions, the
concentrations of metabolites and the actual maximum minimum driving force. The optimization
problem solved is:
```
max min -ΔᵣG
s.t. ΔᵣG = ΔrG⁰ + R T S' ln(C)
     ΔᵣG ≤ 0
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
where `ΔrG` are the Gibbs energies dissipated by the reactions, R is the gas constant, T is
the temperature, S is the stoichiometry of the model, and C is the vector of metabolite
concentrations (and their respective lower and upper bounds).

In case no feasible solution exists, `nothing` is returned.

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
function max_min_driving_force_analysis(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    flux_solution = Dict{String,Float64}(),
    concentration_ratios = Dict{String,Tuple{String,String,Float64}}(),
    proton_ids = ["h_c", "h_e"],
    water_ids = ["h2o_c", "h2o_e"],
    concentration_lb = 1e-9,
    concentration_ub = 100e-3,
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
    modifications = [],
    optimizer,
)
    m = build_max_min_driving_force_model(
        model,
        reaction_standard_gibbs_free_energies;
        flux_solution,
        concentration_lb,
        concentration_ub,
        R,
        T,
        ignore_reaction_ids,
        water_ids,
        proton_ids,
    )
    
    m = add_metabolite_ratio_constraints!(
        m,
        concentration_ratios,
    )
    
    optimized_constraints(
        m;
        objective = m.max_min_driving_force.value,
        optimizer,
        modifications,
    )
end