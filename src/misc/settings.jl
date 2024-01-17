
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

Change the objective sense of optimization. Accepted arguments include
[`Minimal`](@ref), [`Maximal`](@ref), and [`Feasible`](@ref).
"""
set_objective_sense(objective_sense) =
    opt_model -> J.set_objective_sense(opt_model, objective_sense)

export set_objective_sense

"""
$(TYPEDSIGNATURES)

Change the JuMP optimizer used to run the optimization.
"""
set_optimizer(optimizer) = opt_model -> J.set_optimizer(opt_model, optimizer)

export set_optimizer

"""
$(TYPEDSIGNATURES)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
set_optimizer_attribute(attribute_key, value) =
    opt_model -> J.set_optimizer_attribute(opt_model, attribute_key, value)

export set_optimizer_attribute

"""
    silence

Modification that disable all output from the JuMP optimizer (shortcut for
`set_silent` from JuMP).
"""
silence(opt_model) = J.set_silent(opt_model)

"""
$(TYPEDSIGNATURES)

Portable way to set a time limit in seconds for the optimizer computation.
"""
set_time_limit_sec(limit) = opt_model -> J.set_time_limit_sec(opt_model, limit)

export silence
