"""
$(TYPEDSIGNATURES)

Add quasi-thermodynamic constraints to the model to ensure that no
thermodynamically infeasible internal cycles can occur. The function
`metabolite_idx` takes a metabolite ID and the model, and returns the row index
in the stoichiometric matrix (generated using `stoichiometry(model)`)
corresponding to the mass balance around this metabolite. For most models this
will just be `(mid, model) -> first(indexin([mid], metabolites(model)))`. Notable
exceptions include enzyme constrained models, which include virtual enzyme
balances. The latter should not be included in thermodynamic calculations.

Adds the following constraints to the problem:
```
-max_flux_bound × (1 - yᵢ) ≤ xᵢ ≤ max_flux_bound × yᵢ
-max_flux_bound × yᵢ + strict_inequality_tolerance × (1 - yᵢ) ≤ Gᵢ
Gᵢ ≤ -strict_inequality_tolerance × yᵢ + max_flux_bound × (1 - yᵢ)
Nᵢₙₜ' × G = 0
yᵢ ∈ {0, 1}
Gᵢ ∈ ℝ
i ∈ internal reactions
Nᵢₙₜ is the nullspace of the internal stoichiometric matrix (rows =  metabolites, columns = reactions)
```
Note, this modification introduces binary variables, so an optimization solver
capable of handing mixed integer problems needs to be used. The arguments
`max_flux_bound` and `strict_inequality_tolerance` implement the "big-M" method
of indicator constraints.

For more details about the algorithm, see `Schellenberger, Lewis, and, Palsson.
"Elimination of thermodynamically infeasible loops in steady-state metabolic
models.", Biophysical journal, 2011`.
"""
add_loopless_constraints(
    metabolite_idx;
    max_flux_bound = constants.default_reaction_bound, # needs to be an order of magnitude bigger, big M method heuristic
    strict_inequality_tolerance = constants.loopless_strict_inequality_tolerance,
) =
    (model, opt_model) -> begin

        internal_rxn_idxs = [
            ridx for (ridx, rid) in enumerate(reactions(model)) if
            !is_boundary(reaction_stoichiometry(model, rid))
        ]

        metabolite_idxs = metabolite_idx.(metabolites(model), Ref(model)) # TODO need a function that names all the constraints

        N_int =
            nullspace(Array(stoichiometry(model)[metabolite_idxs, internal_rxn_idxs])) # no sparse nullspace function

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

"""
$(TYPEDSIGNATURES)

A convenience wrapper around [`add_loopless_constraints`](@ref), which assumes
that the model being modified only has metabolites in its stoichiometric matrix,
and not, e.g. virtual enzymes.

Supplies `metabolite_idx = (x, m) -> first(indexin([x], metabolites(m)))` to the
base method.
"""
add_loopless_constraints(; kwargs...) =
    add_loopless_constraints((x, m) -> first(indexin([x], metabolites(m))); kwargs...)
