"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the sMOMENT algorithm,
giving a [`SMomentModel`](@ref). The arguments are forwarded to
[`make_smoment_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_smoment(; kwargs...) = model -> make_smoment_model(model; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant which adds extra semantics of the Gecko algorithm,
giving a [`GeckoModel`](@ref). The arguments are forwarded to
[`make_gecko_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_gecko(; kwargs...) = model -> make_gecko_model(model; kwargs...)
