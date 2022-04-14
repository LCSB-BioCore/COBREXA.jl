"""
    protein_dict(model::GeckoModel, opt_model)

Return a dictionary mapping protein molar concentrations to their ids. The
argument `opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref).
"""
protein_dict(model::GeckoModel, opt_model) =
    let gids = genes(model)
        is_solved(opt_model) ?
        Dict(
            [gids[gidx] for (gidx, _) in model.coupling_row_gene_product] .=> _gecko_gene_product_coupling(model) * value.(opt_model[:x]),
        ) : nothing
    end

"""
    protein_dict(model::GeckoModel)

A pipe-able variant of [`protein_dict`](@ref).
"""
protein_dict(model::GeckoModel) = x -> protein_dict(model, x)

"""
    protein_mass_group_dict(model::GeckoModel, opt_model)

Extract the mass utilization in mass groups from a solved [`GeckoModel`](@ref).
"""
protein_mass_group_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(
        (group for (group, _) in model.coupling_row_mass_group) .=>
            _gecko_mass_group_coupling(model) * value.(opt_model[:x]),
    ) : nothing

"""
    protein_mass_group_dict(model::GeckoModel)

A pipe-able variant of [`mass_group_dict`](@ref).
"""
protein_mass_group_dict(model::GeckoModel) = x -> mass_group_dict(model, x)


"""
    protein_mass(model::SMomentModel)

Extract the total mass utilization in a solved [`SMomentModel`](@ref).
"""
protein_mass(model::SMomentModel, opt_model) =
    is_solved(opt_model) ?
    sum((col.capacity_required for col in model.columns) .* value.(opt_model[:x])) : nothing

"""
    protein_mass(model::SMomentModel)


A pipe-able variant of [`protein_mass`](@ref).
"""
protein_mass(model::SMomentModel) = x -> protein_mass(model, x)
