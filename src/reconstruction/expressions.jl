"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the
[`ExpressionLimitedModel`](@ref), simulating the E-Flux algorithm. The
arguments are forwarded to [`make_expression_limited_model`](@ref). Intended
for usage with [`screen`](@ref).
"""
with_expression_limits(; kwargs...) =
    model -> make_expression_limited_model(model; kwargs...)
