"""
$(TYPEDSIGNATURES)

Specifies a model variant with reactions added. Forwards the arguments to
[`add_reactions`](@ref).
"""
with_added_reactions(args...; kwargs...) = m -> add_reactions(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with an added reaction. Forwards the arguments to
[`add_reaction`](@ref).
"""
with_added_reaction(args...; kwargs...) = m -> add_reaction(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with metabolites added. Forwards the arguments to
[`add_metabolites`](@ref).
"""
with_added_metabolites(args...; kwargs...) = m -> add_metabolites(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with an added metabolite. Forwards the arguments to
[`add_metabolite`](@ref).
"""
with_added_metabolite(args...; kwargs...) = m -> add_metabolite(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with genes added. Forwards the arguments to
[`add_genes`](@ref).
"""
with_added_genes(args...; kwargs...) = m -> add_genes(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with an added gene. Forwards the arguments to
[`add_gene`](@ref).
"""
with_added_gene(args...; kwargs...) = m -> add_gene(m, args...; kwargs...)

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
[`remove_metabolite`](@ref).
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

Specifies a model variant without a certain reaction. Forwards arguments to
[`remove_reaction`](@ref).
"""
with_removed_reaction(args...; kwargs...) = m -> remove_reaction(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Plural version of [`with_removed_reaction`](@ref), calls
[`remove_reactions`](@ref) internally.
"""
with_removed_reactions(args...; kwargs...) = m -> remove_reactions(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Species a model variant that adds a biomass metabolite to the biomass reaction.
Forwards arguments to [`add_biomass_metabolite`](@ref).
"""
with_added_biomass_metabolite(args...; kwargs...) =
    m -> add_biomass_metabolite(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Specifies a model variant with the bounds changed for the gene product. Forwards
the arguments to [`change_gene_product_bound`](@ref).
"""
with_changed_gene_product_bound(args...; kwargs...) =
    m -> change_gene_product_bound(m, args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Plural version of [`with_changed_gene_product_bound`](@ref), calls
[`change_gene_product_bounds`](@ref) internally.
"""
with_changed_gene_product_bounds(args...; kwargs...) =
    m -> change_gene_product_bounds(m, args...; kwargs...)
