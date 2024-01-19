
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

Make an JuMP model out of `constraints` using [`optimization_model`](@ref)
(most arguments are forwarded there), then apply the `settings`, optimize
the model, and return either `nothing` if the optimization failed, or `output`
substituted with the solved values (`output` defaults to `constraints`.

For a "nice" version for simpler finding of metabolic model optima, use
[`flux_balance`](@ref).
"""
function optimized_constraints(
    constraints::C.ConstraintTree;
    settings = [],
    output::C.ConstraintTreeElem = constraints,
    kwargs...,
)
    om = optimization_model(constraints; kwargs...)
    for m in settings
        m(om)
    end

    optimized_model(om; output)
end

export optimized_constraints
