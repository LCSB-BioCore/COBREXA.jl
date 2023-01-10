"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the sMOMENT algorithm,
giving a [`SimplifiedEnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_simplified_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_simplified_enzyme_constrained(; kwargs...) =
    model -> make_simplified_enzyme_constrained_model(model; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the EnzymeConstrained algorithm,
giving a [`EnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_enzyme_constrained(; kwargs...) =
    model -> make_enzyme_constrained_model(model; kwargs...)
