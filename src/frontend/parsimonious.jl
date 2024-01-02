
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

Optimize the system of `constraints` to get the optimal `objective` value. Then
try to find a "parsimonious" solution with the same `objective` value, which
optimizes the `parsimonious_objective` (possibly also switching optimization
sense, optimizer, and adding more modifications).

For efficiency, everything is performed on a single instance of JuMP model.

A simpler version suitable for direct work with metabolic models is available
in [`parsimonious_flux_balance`](@ref).
"""
function parsimonious_optimized_constraints(
    constraints::C.ConstraintTreeElem;
    objective::C.Value,
    modifications = [],
    parsimonious_objective::C.Value,
    parsimonious_optimizer = nothing,
    parsimonious_sense = J.MIN_SENSE,
    parsimonious_modifications = [],
    tolerances = [absolute_tolerance_bound(0)],
    output = constraints,
    kwargs...,
)

    # first solve the optimization problem with the original objective
    om = optimization_model(constraints; objective, kwargs...)
    for m in modifications
        m(om)
    end
    J.optimize!(om)
    is_solved(om) || return nothing

    target_objective_value = J.objective_value(om)

    # switch to parsimonizing the solution w.r.t. to the objective value
    isnothing(parsimonious_optimizer) || J.set_optimizer(om, parsimonious_optimizer)
    for m in parsimonious_modifications
        m(om)
    end

    J.@objective(om, J.MIN_SENSE, C.substitute(parsimonious_objective, om[:x]))

    # try all admissible tolerances
    for tolerance in tolerances
        (lb, ub) = tolerance(target_objective_value)
        J.@constraint(
            om,
            pfba_tolerance_constraint,
            lb <= C.substitute(objective, om[:x]) <= ub
        )

        J.optimize!(om)
        is_solved(om) && return C.constraint_values(output, J.value.(om[:x]))

        J.delete(om, pfba_tolerance_constraint)
        J.unregister(om, :pfba_tolerance_constraint)
    end

    # all tolerances failed
    return nothing
end

export parsimonious_optimized_constraints

"""
$(TYPEDSIGNATURES)

Compute a parsimonious flux solution for the given `model`. In short, the
objective value of the parsimonious solution should be the same as the one from
[`flux_balance_analysis`](@ref), except the squared sum of reaction fluxes is minimized.
If there are multiple possible fluxes that achieve a given objective value,
parsimonious flux thus represents the "minimum energy" one, thus arguably more
realistic. The optimized squared distance is present in the result as
`parsimonious_objective`.

Most arguments are forwarded to [`parsimonious_optimized_constraints`](@ref),
with some (objectives) filled in automatically to fit the common processing of
FBC models, and some (`tolerances`) provided with more practical defaults.

Similarly to the [`flux_balance_analysis`](@ref), returns a tree with the optimization
solutions of the shape as given by [`fbc_model_constraints`](@ref).
"""
function parsimonious_flux_balance_analysis(
    model::A.AbstractFBCModel,
    optimizer;
    tolerances = relative_tolerance_bound.(1 .- [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]),
    kwargs...,
)
    constraints = fbc_model_constraints(model)
    parsimonious_objective = squared_sum_objective(constraints.fluxes)
    parsimonious_optimized_constraints(
        constraints * :parsimonious_objective^C.Constraint(parsimonious_objective);
        optimizer,
        objective = constraints.objective.value,
        parsimonious_objective,
        tolerances,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Pipe-able variant of [`parsimonious_flux_balance_analysis`](@ref).
"""
parsimonious_flux_balance_analysis(optimizer; kwargs...) =
    model -> parsimonious_flux_balance_analysis(model, optimizer; kwargs...)

export parsimonious_flux_balance_analysis
