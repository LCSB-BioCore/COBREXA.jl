
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

Optimize the system given by `constraints` and `objective` with `optimizer`
(with custom `settings`) for all combination of constriants given by `dims`.

`dims` should be compatible with pairs that assign a sequence of breaks to a
`ConstraintTrees.Value`: For example, `organism.fluxes.PFK => 1:3` will compute
optima of the model with the flux through `PFK` constrained to be equal to 1, 2
and 3.

In turn, all `dims` are converted to groups of equality constraints, and the
model is solved for all combinations. Shape of the output matrix corresponds to
`Iterators.product(last.(dims)...)`.

Operation is parallelized by distribution over `workers`; by default all
`Distributed` workers are used.
"""
function constraints_objective_envelope(
    constraints::C.ConstraintTree,
    dims...;
    objective::C.Value,
    sense = Maximal,
    optimizer,
    settings = [],
    workers = D.workers(),
)
    values = first.(dims)
    ranges = last.(dims)

    screen_optimization_model(
        constraints,
        Iterators.product(ranges...);
        objective,
        sense,
        optimizer,
        settings,
        workers,
    ) do om, coords
        con_refs = [
            begin
                J.@constraint(om, con_ref, C.substitute(v, om[:x]) == x)
                con_ref
            end for (v, x) in zip(values, coords)
        ]
        J.optimize!(om)
        res = is_solved(om) ? J.objective_value(om) : nothing
        J.delete.(con_refs)
        res
    end
end
