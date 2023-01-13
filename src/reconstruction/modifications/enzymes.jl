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

"""
$(TYPEDSIGNATURES)

Specifies a model variant that adds a pseudoribosome to a model. Args and kwargs
are forwarded to [`add_pseudoribosome`](@ref).
"""
with_pseudoribosome(args...; kwargs...) =
    model -> add_pseudoribosome(model, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant that overwrites the current isozymes associated with the
model through calling [`add_isozymes`](@ref).
"""
with_isozymes(args...) = model -> add_isozymes(model, args...)
