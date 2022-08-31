"""
$(TYPEDSIGNATURES)

Return a dictionary mapping protein molar concentrations to their ids. The
argument `opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref). See [`flux_dict`](@ref) for the corresponding
function that returns a dictionary of solved fluxes.
"""
gene_product_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(genes(model) .=> value.(opt_model[:x])[(length(model.columns)+1):end]) : nothing

"""
$(TYPEDSIGNATURES)

A pipe-able variant of [`gene_product_dict`](@ref).
"""
gene_product_dict(model::GeckoModel) = x -> gene_product_dict(model, x)

"""
$(TYPEDSIGNATURES)

Extract the mass utilization in mass groups from a solved [`GeckoModel`](@ref).
"""
gene_product_mass_group_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(
        grp.group_id => dot(
            value.(opt_model[:x])[length(model.columns).+grp.gene_product_idxs],
            grp.gene_product_molar_masses,
        ) for grp in model.coupling_row_mass_group
    ) : nothing

"""
$(TYPEDSIGNATURES)

A pipe-able variant of [`gene_product_mass_group_dict`](@ref).
"""
gene_product_mass_group_dict(model::GeckoModel) =
    x -> gene_product_mass_group_dict(model, x)

"""
$(TYPEDSIGNATURES)

Extract the total mass utilization in a solved [`SMomentModel`](@ref).
"""
gene_product_mass(model::SMomentModel, opt_model) =
    is_solved(opt_model) ?
    sum((col.capacity_required for col in model.columns) .* value.(opt_model[:x])) : nothing

"""
$(TYPEDSIGNATURES)

A pipe-able variant of [`gene_product_mass`](@ref).
"""
gene_product_mass(model::SMomentModel) = x -> gene_product_mass(model, x)
