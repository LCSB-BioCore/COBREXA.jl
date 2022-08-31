"""
$(TYPEDSIGNATURES)

Specifies a model variant that has a new bound set. Forwards arguments to
[`change_bound`](@ref). Intended for usage with [`screen`](@ref).
"""
with_changed_bound(args...; kwargs...) = m -> change_bound(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant that has new bounds set. Forwards arguments to
[`change_bounds`](@ref). Intended for usage with [`screen`](@ref).
"""
with_changed_bounds(args...; kwargs...) = m -> change_bounds(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant without a certain metabolite. Forwards arguments to
[`remove_metabolite`](@ref). Intended to be used with [`screen`](@ref).
"""
with_removed_metabolite(args...; kwargs...) = m -> remove_metabolite(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Plural version of [`with_removed_metabolite`](@ref), calls
[`remove_metabolites`](@ref) internally.
"""
with_removed_metabolites(args...; kwargs...) =
    m -> remove_metabolites(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with reactions added. Forwards the arguments to
[`add_reactions`](@ref). Intended to be used with [`screen`](@ref).
"""
with_added_reactions(args...; kwargs...) = m -> add_reactions(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant without a certain reaction. Forwards arguments to
[`remove_reaction`](@ref). Intended to be used with [`screen`](@ref).
"""
with_removed_reaction(args...; kwargs...) = m -> remove_reaction(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Plural version of [`with_removed_reaction`](@ref), calls
[`remove_reactions`](@ref) internally.
"""
with_removed_reactions(args...; kwargs...) = m -> remove_reactions(m, args...; kwargs...)
