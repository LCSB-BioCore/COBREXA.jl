
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

Linearly scale all bounds in a constraint tree by the `factor`. This actually
changes the model semantics, and may not work in surprising/improper ways with
some constraint systems, esp. the MILP and QP ones.

See also [`scale_constraints`](@ref).
"""
scale_bounds(tree::C.ConstraintTree, factor) =
    C.map(tree) do c
        isnothing(c.bound) ? c : C.Constraint(value = c.value, bound = factor * c.bound)
    end

"""
$(TYPEDSIGNATURES)

Linearly scale all constraints in a constraint tree by the `factor`.

See also [`scale_bounds`](@ref).
"""
scale_constraints(tree::C.ConstraintTree, factor) =
    C.map(tree) do c
        c * factor
    end
