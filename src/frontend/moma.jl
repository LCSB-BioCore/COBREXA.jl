
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

Find a feasible solution of the "minimal metabolic adjustment analysis" (MOMA)
for the `model`, which is the "closest" feasible solution to the given
`reference_fluxes`, in the sense of squared-sum error distance. The minimized
squared distance (the objective) is present in the result tree as
`minimal_adjustment_objective`.

This is often used for models with smaller feasible region than the reference
models (typically handicapped by a knockout, nutritional deficiency or a
similar perturbation). MOMA solution then gives an expectable "easiest"
adjustment of the organism towards a somewhat working state.

Reference fluxes that do not exist in the model are ignored (internally, the
objective is constructed via [`squared_sum_error_value`](@ref)).

Additional parameters are forwarded to [`optimized_constraints`](@ref).
"""
function minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_fluxes::Dict{Symbol,Float64},
    optimizer;
    kwargs...,
)
    constraints = fbc_flux_balance_constraints(model)
    objective = squared_sum_error_value(constraints.fluxes, reference_fluxes)
    optimized_constraints(
        constraints * :minimal_adjustment_objective^C.Constraint(objective);
        optimizer,
        objective,
        sense = Minimal,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

A slightly easier-to-use version of
[`minimization_of_metabolic_adjustment_analysis`](@ref) that computes the
reference flux as the optimal solution of the [`reference_model`](@ref). The
reference flux is calculated using `reference_optimizer` and
`reference_modifications`, which default to the `optimizer` and `settings`.

Leftover arguments are passed to the overload of
[`minimization_of_metabolic_adjustment_analysis`](@ref) that accepts the
reference flux dictionary.
"""
function minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel,
    optimizer;
    reference_optimizer = optimizer,
    settings = [],
    reference_settings = settings,
    kwargs...,
)
    reference_constraints = fbc_flux_balance_constraints(reference_model)
    reference_fluxes = optimized_constraints(
        reference_constraints;
        optimizer = reference_optimizer,
        settings = reference_settings,
        output = reference_constraints.fluxes,
    )
    isnothing(reference_fluxes) && return nothing
    minimization_of_metabolic_adjustment_analysis(
        model,
        reference_fluxes,
        optimizer;
        settings,
        kwargs...,
    )
end

export minimization_of_metabolic_adjustment_analysis

"""
$(TYPEDSIGNATURES)

Like [`minimization_of_metabolic_adjustment_analysis`](@ref) but optimizes the
L1 norm. This typically produces a sufficiently good result with less
resources, depending on the situation. See documentation of
[`linear_parsimonious_flux_balance_analysis`](@ref) for some of the
considerations.
"""
function linear_minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_fluxes::Dict{Symbol,Float64},
    optimizer;
    kwargs...,
)
    constraints = fbc_flux_balance_constraints(model)

    difference = C.zip(ct.fluxes, C.Tree(reference_fluxes)) do orig, ref
        C.Constraint(orig.value - ref)
    end

    difference_split_variables =
        C.variables(keys = keys(difference), bounds = C.Between(0, Inf))
    constraints += :reference_positive_diff^difference_split_variables
    constraints += :reference_negative_diff^difference_split_variables

    # `difference` actually doesn't need to go to the CT, but we include it
    # anyway to calm the curiosity of good neighbors.
    constraints *= :reference_diff^difference
    constraints *=
        :reference_directional_diff_balance^sign_split_constraints(
            positive = constraints.reference_positive_diff,
            negative = constraints.reference_negative_diff,
            signed = difference,
        )

    objective =
        sum_value(constraints.reference_positive_diff, constraints.reference_negative_diff)

    optimized_constraints(
        constraints * :linear_minimal_adjustment_objective^C.Constraint(objective);
        optimizer,
        objective,
        sense = Minimal,
        kwargs...,
    )
end

function linear_minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel,
    optimizer;
    reference_optimizer = optimizer,
    settings = [],
    reference_settings = settings,
    kwargs...,
)
    reference_constraints = fbc_flux_balance_constraints(reference_model)
    reference_fluxes = optimized_constraints(
        reference_constraints;
        optimizer = reference_optimizer,
        settings = reference_settings,
        output = reference_constraints.fluxes,
    )
    isnothing(reference_fluxes) && return nothing
    linear_minimization_of_metabolic_adjustment_analysis(
        model,
        reference_fluxes,
        optimizer;
        settings,
        kwargs...,
    )
end

export linear_minimization_of_metabolic_adjustment_analysis
