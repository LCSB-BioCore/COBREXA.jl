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
