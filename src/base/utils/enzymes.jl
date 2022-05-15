"""
    gene_product_dict(model::GeckoModel, opt_model)

Return a dictionary mapping protein molar concentrations to their ids. The
argument `opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref). See [`flux_dict`](@ref) for the corresponding
function that returns a dictionary of solved fluxes.
"""
gene_product_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(genes(model) .=> value.(opt_model[:x])[(end-n_genes(model)+1):end]) : nothing

"""
    gene_product_dict(model::GeckoModel)

A pipe-able variant of [`gene_product_dict`](@ref).
"""
gene_product_dict(model::GeckoModel) = x -> gene_product_dict(model, x)

"""
    gene_product_mass_group_dict(model::GeckoModel, opt_model)

Extract the mass utilization in mass groups from a solved [`GeckoModel`](@ref).
"""
gene_product_mass_group_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(
        grp.group_id => dot(
            value.(opt_model[:x])[(n_reactions(model)-n_genes(model)).+grp.gene_product_idxs],
            grp.gene_product_molar_masses,
        ) for grp in model.coupling_row_mass_group
    ) : nothing

"""
    gene_product_mass_group_dict(model::GeckoModel)

A pipe-able variant of [`mass_group_dict`](@ref).
"""
gene_product_mass_group_dict(model::GeckoModel) = x -> mass_group_dict(model, x)

"""
    gene_product_mass(model::SMomentModel)

Extract the total mass utilization in a solved [`SMomentModel`](@ref).
"""
gene_product_mass(model::SMomentModel, opt_model) =
    is_solved(opt_model) ?
    sum((col.capacity_required for col in model.columns) .* value.(opt_model[:x])) : nothing

"""
    gene_product_mass(model::SMomentModel)


A pipe-able variant of [`gene_product_mass`](@ref).
"""
gene_product_mass(model::SMomentModel) = x -> gene_product_mass(model, x)
