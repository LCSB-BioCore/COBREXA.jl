"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the sMOMENT algorithm,
giving a [`SimplifiedEnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_simplified_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_simplified_enzyme_constraints(; kwargs...) =
    model -> make_simplified_enzyme_constrained_model(model; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the EnzymeConstrained algorithm,
giving a [`EnzymeConstrainedModel`](@ref). The arguments are forwarded to
[`make_enzyme_constrained_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_enzyme_constraints(; kwargs...) =
    model -> make_enzyme_constrained_model(model; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant that adds a virtualribosome to a model. Args and kwargs
are forwarded to [`add_virtualribosome`](@ref).
"""
with_virtualribosome(args...; kwargs...) =
    model -> add_virtualribosome(model, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant that overwrites the current isozymes associated with the
model through calling [`add_isozymes`](@ref).
"""
with_isozymes(args...; kwargs...) = model -> add_isozymes(model, args...; kwargs...)
