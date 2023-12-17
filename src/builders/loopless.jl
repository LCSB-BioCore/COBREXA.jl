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
add_loopless_constraints(;
    max_flux_bound = _constants.default_reaction_bound, # needs to be an order of magnitude bigger, big M method heuristic
    strict_inequality_tolerance = _constants.loopless_strict_inequality_tolerance,
) =
    (model, opt_model) -> begin

        internal_rxn_idxs = [
            ridx for (ridx, rid) in enumerate(reactions(model)) if
            !is_boundary(reaction_stoichiometry(model, rid))
        ]

        N_int = nullspace(Array(stoichiometry(model)[:, internal_rxn_idxs])) # no sparse nullspace function

        y = @variable(opt_model, y[1:length(internal_rxn_idxs)], Bin)
        G = @variable(opt_model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

        x = opt_model[:x]
        for (cidx, ridx) in enumerate(internal_rxn_idxs)
            @constraint(opt_model, -max_flux_bound * (1 - y[cidx]) <= x[ridx])
            @constraint(opt_model, x[ridx] <= max_flux_bound * y[cidx])

            @constraint(
                opt_model,
                -max_flux_bound * y[cidx] + strict_inequality_tolerance * (1 - y[cidx]) <= G[cidx]
            )
            @constraint(
                opt_model,
                G[cidx] <=
                -strict_inequality_tolerance * y[cidx] + max_flux_bound * (1 - y[cidx])
            )
        end

        @constraint(opt_model, N_int' * G .== 0)
    end
