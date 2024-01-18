
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

TODO
"""
function constraints_variability(
    constraints::C.ConstraintTree,
    targets::C.ConstraintTree;
    optimizer,
    settings = [],
    workers = D.workers(),
)::C.Tree{Tuple{Maybe{Float64},Maybe{Float64}}}

    #TODO settings?
    target_array = [dim * dir for dim in tree_deflate(C.value, x), dir in (-1, 1)]

    result_array = screen_optimization_model(
        constraints,
        target_array;
        optimizer,
        settings,
        workers,
    ) do om, target
        J.@objective(om, Maximal, C.substitute(target, om[:x]))
        J.optimize!(om)
        is_solved(om) ? J.objective_value(om) : nothing
    end

    constraint_tree_reinflate(targets, [tuple(a, b) for (a, b) in eachrow(result_array)])
end
