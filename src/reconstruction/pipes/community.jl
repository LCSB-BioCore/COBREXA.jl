"""
$(TYPEDSIGNATURES)

Specifies a model variant where the abundances of community members has been
changed. Forwards arguments to [`change_abundances`](@ref).
"""
with_changed_abundances(args...; kwargs...) = m -> change_abundances(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant where an environmental exchange reaction has its
bounds changed. Calls [`change_environmental_bound`](@ref) internally.
"""
with_changed_environmental_bound(args...; kwargs...) =
    m -> change_environmental_bound(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Plural variant of [`with_changed_environmental_bound`](@ref)
"""
with_changed_environmental_bounds(args...; kwargs...) =
    m -> change_environmental_bounds(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Species a model variant that wraps a [`CommunityModel`](@ref) into a
[`EqualGrowthCommunityModel`](@ref). Forwards the arguments to the constructor.
"""
with_equal_growth_objective(args...; kwargs...) = m -> make_EqualGrowthCommunityModel(m, args...; kwargs...)