"""
    protein_dict(model::GeckoModel, opt_model)

Return a dictionary mapping protein molar concentrations to their ids. The
argument `opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref). See [`flux_dict`](@ref) for the corresponding
function that returns a dictionary of solved fluxes.
"""
protein_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(genes(model) .=> value.(opt_model[:x])[(n_reactions(model)+1):end]) : nothing

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
        grp[1] => dot(value.(opt_model[:x])[n_reactions(model).+grp[2]], grp[3]) for
        grp in model.coupling_row_mass_group
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
