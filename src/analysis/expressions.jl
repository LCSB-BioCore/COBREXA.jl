
"""
$(TYPEDSIGNATURES)

Build an [`ExpressionLimitedModel`](@ref).
"""
make_expression_limited_model(
    model::MetabolicModel;
    relative_expression::Dict{String,Float64},
    bounding_function::Function = expression_probabilistic_bounds,
) = ExpressionLimitedModel(; relative_expression, bounding_function, inner = model)
