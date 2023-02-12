"""
$(TYPEDSIGNATURES)

Extract the mass utilization in mass groups from a solved [`EnzymeConstrainedModel`](@ref).
"""
gene_product_mass_group_dict(model::EnzymeConstrainedModel, opt_model) =
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
gene_product_mass_group_dict(model::EnzymeConstrainedModel) =
    x -> gene_product_mass_group_dict(model, x)

"""
$(TYPEDSIGNATURES)

Extract the total mass utilization in a solved [`SimplifiedEnzymeConstrainedModel`](@ref).
"""
gene_product_mass(model::SimplifiedEnzymeConstrainedModel, opt_model) =
    is_solved(opt_model) ?
    sum((col.capacity_required for col in model.columns) .* value.(opt_model[:x])) : nothing

"""
$(TYPEDSIGNATURES)

A pipe-able variant of [`gene_product_mass`](@ref).
"""
gene_product_mass(model::SimplifiedEnzymeConstrainedModel) =
    x -> gene_product_mass(model, x)
