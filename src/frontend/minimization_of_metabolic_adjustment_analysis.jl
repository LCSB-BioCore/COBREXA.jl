
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
objective is constructed via [`squared_sum_error_objective`](@ref)).

Additional parameters are forwarded to [`optimized_constraints`](@ref).
"""
function minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_fluxes::Dict{Symbol,Float64},
    optimizer;
    kwargs...,
)
    constraints = fbc_model_constraints(model)
    objective = squared_sum_error_objective(constraints.fluxes, reference_fluxes)
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

A slightly easier-to-use version of [`minimization_of_metabolic_adjustment_analysis`](@ref) that
computes the reference flux as the optimal solution of the
[`reference_model`](@ref). The reference flux is calculated using
`reference_optimizer` and `reference_modifications`, which default to the
`optimizer` and `modifications`.

Leftover arguments are passed to the overload of
[`minimization_of_metabolic_adjustment_analysis`](@ref) that accepts the reference flux
dictionary.
"""
function minimization_of_metabolic_adjustment_analysis(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel,
    optimizer;
    reference_optimizer = optimizer,
    modifications = [],
    reference_modifications = modifications,
    kwargs...,
)
    reference_constraints = fbc_model_constraints(reference_model)
    reference_fluxes = optimized_constraints(
        reference_constraints;
        optimizer = reference_optimizer,
        modifications = reference_modifications,
        output = reference_constraints.fluxes,
    )
    isnothing(reference_fluxes) && return nothing
    minimization_of_metabolic_adjustment_analysis(
        model,
        reference_fluxes,
        optimizer;
        modifications,
        kwargs...,
    )
end

export minimization_of_metabolic_adjustment_analysis
