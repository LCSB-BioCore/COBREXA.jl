# constructors for enzyme constrained models

"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the sMOMENT algorithm,
giving a [`SimplifiedEnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_simplified_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_simplified_enzyme_constraints(args...; kwargs...) =
    model -> make_simplified_enzyme_constrained_model(model, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the EnzymeConstrained algorithm,
giving a [`EnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_enzyme_constraints(args...; kwargs...) =
    model -> make_enzyme_constrained_model(model, args...; kwargs...)
