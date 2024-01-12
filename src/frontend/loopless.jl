
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

Perform a flux balance analysis with added quasi-thermodynamic constraints that
ensure that thermodynamically infeasible internal cycles can not occur. The
method is closer described by: Schellenberger, Lewis, and, Palsson.
"Elimination of thermodynamically infeasible loops in steady-state metabolic
models.", Biophysical journal, 2011`.

The loopless condition comes with a performance overhead: the computation needs
to find the null space of the stoichiometry matrix (essentially inverting it);
and the subsequently created optimization problem contains binary variables for
each internal reaction, thus requiring a MILP solver and a potentially
exponential solving time.

The arguments `driving_force_max_bound` and `driving_force_nonzero_bound` set
the bounds (possibly negated ones) on the virtual "driving forces" (G_i in the
paper).
"""
function loopless_flux_balance_analysis(
    model::A.AbstractFBCModel;
    flux_infinity_bound = 10000.0,
    driving_force_nonzero_bound = 1.0,
    driving_force_infinity_bound = 1000.0,
    settings = [],
    optimizer,
)

    constraints = flux_balance_constraints(model)

    rxns = A.reactions(model)
    stoi = A.stoichiometry(model)
    internal_mask = count(stoi .!= 0; dims = 1)[begin, :] .> 1
    internal_reactions =
        Symbol.(rxns[reactions_internal]), constraints =
            constraints +
            :loopless_directions^C.variables(
                keys = internal_reactions,
                bounds = Switch(0, 1),
            ) +
            :loopless_driving_forces^C.variables(keys = internal_reactions)

    constraints *=
        :loopless_constraints^loopless_constraints(;
            fluxes = constraints.fluxes,
            loopless_directions = constraints.loopless_directions,
            loopless_driving_forces = constraints.loopless_driving_forces,
            internal_reactions,
            internal_nullspace = LinearAlgebra.nullspace(Matrix(stoi[:, internal_mask])),
            flux_infinity_bound,
            driving_force_nonzero_bound,
            driving_force_infinity_bound,
        )

    optimized_constraints(m; objective = m.objective.value, optimizer, settings)
end

export loopless_flux_balance_analysis
