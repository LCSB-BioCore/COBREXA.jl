"""
    with_set_bound(args...; kwargs...) 

Specifies a model variant that has a new bound set. Forwards arguments to
[`set_bound`](@ref). Intended for usage with [`screen`](@ref).
"""
with_set_bound(args...; kwargs...) = m -> set_bound(m, args...; kwargs...)

"""
    with_removed_metabolites(args...; kwargs...) 

Specifies a model variant without specified metabolites. Forwards arguments to
[`remove_metabolites`](@ref). Intended to be used with [`screen`](@ref).
"""
with_removed_metabolites(args...; kwargs...) =
    m -> remove_metabolites(m, args...; kwargs...)


"""
    with_added_reactions(args...; kwargs...) 

Specifies a model variant with reactions added. Forwards the arguments to
[`add_reactions`](@ref). Intended to be used with [`screen`](@ref).
"""
with_added_reactions(args...; kwargs...) = m -> add_reactions(m, args...; kwargs...)

"""
    with_removed_reactions(args...; kwargs...) 

Specifies a model variant with specified reactions removed. Forwards arguments
to [`remove_reactions`](@ref). Intended to be used with [`screen`](@ref).
"""
with_removed_reactions(args...; kwargs...) = m -> remove_reactions(m, args...; kwargs...)
