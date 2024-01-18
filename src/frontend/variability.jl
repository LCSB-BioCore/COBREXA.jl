
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
function flux_variability_analysis(
    model::A.AbstractFBCModel;
    objective_bound,
    optimizer,
    settings,
    workers = D.workers(),
)
    constraints = flux_balance_constraints(model)

    objective = constraints.objective_value

    objective_flux = optimized_constraints(
        constraints;
        objective = constraints.objective.value,
        output = constraints.objective,
        optimizer,
        settings,
    )

    isnothing(objective_flux) && return nothing

    constraint_variability(
        constraints *
        :objective_bound^C.Constraint(objective, objective_bound(objective_flux)),
        constraints.fluxes;
        optimizer,
        settings,
        workers,
    )
end

export flux_variability_analysis
