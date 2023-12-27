
"""
$(TYPEDSIGNATURES)

Add loopless constraints to the model, `m`. Specify the internal reactions with
`internal_reaction_ids`, as well as the columns of the stoichiometric nullspace
of these reactions in `internal_reaction_stoichiometry_nullspace_columns`. See
the example below for more information about the nullspace. By default, the
`fluxes` of the model `m` are made loopless.

The Big-M method is used to ensure the sign of fluxes and pseudo Gibbs free
energy of reactions match. For this, ensure that `max_flux_bound` is at least
one order of magnitude bigger than the largest expected absolute flux value.
Additionally, ensure that `strict_inequality_tolerance` is smaller than any
expected pseudo Gibbs free energy of reaction value. The defaults work well for
most problems.

# Example
```
internal_reaction_stoichiometry_nullspace_columns =
    eachcol(nullspace(Array(A.stoichiometry(model)[:, internal_rxn_idxs_in_order_of_internal_rxn_ids])))
```
"""
function add_loopless_constraints!(
    m,
    internal_reaction_ids,
    internal_reaction_stoichiometry_nullspace_columns;
    fluxes = m.fluxes,
    max_flux_bound = 1000.0, # needs to be an order of magnitude bigger, big M method heuristic
    strict_inequality_tolerance = 1.0, # heuristic from paper
)

    # add loopless variables: flux direction setters (binary) and pseudo gibbs free energy of reaction
    m +=
        :loopless_binary_variables^C.variables(
            keys = internal_reaction_ids,
            bounds = Binary(),
        )
    m +=
        :pseudo_gibbs_free_energy_reaction^C.variables(
            keys = internal_reaction_ids,
            bounds = C.Between(-Inf, Inf),
        )

    # add -1000 * (1-a) ≤ v ≤ 1000 * a which need to be split into forward and backward components
    # -1000 * (1-a) - v ≤ 0 (backward)
    # v - 1000a ≤ 0 (forward)
    m *=
        :loopless_reaction_directions^:backward^C.ConstraintTree(
            rid => C.Constraint(
                value = -max_flux_bound * (1 - m.loopless_binary_variables[rid].value) -
                        fluxes[rid].value,
                bound = C.Between(Inf, 0),
            ) for rid in internal_reaction_ids
        )
    m *=
        :loopless_reaction_directions^:forward^C.ConstraintTree(
            rid => C.Constraint(
                value = fluxes[rid].value -
                        max_flux_bound * m.loopless_binary_variables[rid].value,
                bound = C.Between(-Inf, 0),
            ) for rid in internal_reaction_ids
        )

    # add -1000*a + 1 * (1-a) ≤ Gibbs ≤ -1 * a + 1000 * (1 - a) which also need to be split
    # -1000 * a + 1 * (1-a)  - G ≤ 0 (backward)
    # G + 1 * a - 1000 * (1-a) ≤ 0 (forward)
    m *=
        :loopless_pseudo_gibbs_sign^:backward^C.ConstraintTree(
            rid => C.Constraint(
                value = -max_flux_bound * m.loopless_binary_variables[rid].value +
                        strict_inequality_tolerance *
                        (1 - m.loopless_binary_variables[rid].value) -
                        m.pseudo_gibbs_free_energy_reaction[rid].value,
                bound = (-Inf, 0),
            ) for rid in internal_reaction_ids
        )
    m *=
        :loopless_pseudo_gibbs_sign^:forward^C.ConstraintTree(
            rid => C.Constraint(
                value = m.pseudo_gibbs_free_energy_reaction[rid].value +
                        strict_inequality_tolerance *
                        m.loopless_binary_variables[rid].value -
                        max_flux_bound * (1 - m.loopless_binary_variables[rid].value),
                bound = (-Inf, 0),
            ) for rid in internal_reaction_ids
        )

    # use nullspace to ensure there are no loops 
    m *=
        :loopless_condition^C.ConstraintTree(
            Symbol(:nullspace_vector, i) => C.Constraint(
                value = sum(
                    col[j] *
                    m.pseudo_gibbs_free_energy_reaction[internal_reaction_ids[j]].value
                    for j in eachindex(col)
                ),
                bound = C.EqualTo(0),
            ) for (i, col) in enumerate(internal_reaction_stoichiometry_nullspace_columns)
        )

    m
end

export add_loopless_constraints!

"""
$(TYPEDSIGNATURES)

Add quasi-thermodynamic constraints to the model to ensure that no thermodynamically
infeasible internal cycles can occur. Adds the following constraints to the problem:
```
-max_flux_bound × (1 - yᵢ) ≤ xᵢ ≤ max_flux_bound × yᵢ
-max_flux_bound × yᵢ + strict_inequality_tolerance × (1 - yᵢ) ≤ Gᵢ
Gᵢ ≤ -strict_inequality_tolerance × yᵢ + max_flux_bound × (1 - yᵢ)
Nᵢₙₜ' × G = 0
yᵢ ∈ {0, 1}
Gᵢ ∈ ℝ
i ∈ internal reactions
Nᵢₙₜ is the nullspace of the internal stoichiometric matrix
```
Note, this modification introduces binary variables, so an optimization solver capable of
handing mixed integer problems needs to be used. The arguments `max_flux_bound` and
`strict_inequality_tolerance` implement the "big-M" method of indicator constraints.

For more details about the algorithm, see `Schellenberger, Lewis, and, Palsson. "Elimination
of thermodynamically infeasible loops in steady-state metabolic models.", Biophysical
journal, 2011`.
"""
function loopless_flux_balance_analysis(
    model;
    max_flux_bound = 1000.0, # needs to be an order of magnitude bigger, big M method heuristic
    strict_inequality_tolerance = 1.0, # heuristic from paper
    modifications = [],
    optimizer,
)

    m = fbc_model_constraints(model)

    # find all internal reactions
    internal_reactions = [
        (i, Symbol(rid)) for
        (i, rid) in enumerate(A.reactions(model)) if !is_boundary(model, rid)
    ]
    internal_reaction_ids = last.(internal_reactions)
    internal_reaction_idxs = first.(internal_reactions) # order needs to match the internal reaction ids below

    internal_reaction_stoichiometry_nullspace_columns =
        eachcol(nullspace(Array(A.stoichiometry(model)[:, internal_reaction_idxs]))) # no sparse nullspace function

    m = add_loopless_constraints!(
        m,
        internal_reaction_ids,
        internal_reaction_stoichiometry_nullspace_columns;
        max_flux_bound,
        strict_inequality_tolerance,
    )

    # solve
    optimized_constraints(m; objective = m.objective.value, optimizer, modifications)
end

export loopless_flux_balance_analysis
