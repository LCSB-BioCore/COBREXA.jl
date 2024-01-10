
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
    settings = [],
    optimizer,
)

    m = fbc_flux_balance_constraints(model)

    # find all internal reactions
    internal_reactions = [
        (i, Symbol(rid)) for
        (i, rid) in enumerate(A.reactions(model)) if !is_boundary(model, rid)
    ]
    internal_reaction_ids = last.(internal_reactions)
    internal_reaction_idxs = first.(internal_reactions) # order needs to match the internal reaction ids below

    internal_reaction_stoichiometry_nullspace_columns = eachcol(
        LinearAlgebra.nullspace(Array(A.stoichiometry(model)[:, internal_reaction_idxs])),
    ) # no sparse nullspace function

    m = with_loopless_constraints(
        m,
        internal_reaction_ids,
        internal_reaction_stoichiometry_nullspace_columns;
        max_flux_bound,
        strict_inequality_tolerance,
    )

    # solve
    optimized_constraints(m; objective = m.objective.value, optimizer, settings)
end

export loopless_flux_balance_analysis
