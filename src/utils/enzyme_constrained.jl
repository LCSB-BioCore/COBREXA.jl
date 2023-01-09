"""
$(TYPEDSIGNATURES)

Return a dictionary mapping protein molar concentrations to their ids. The
argument `opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref). See [`flux_dict`](@ref) for the corresponding
function that returns a dictionary of solved fluxes.
"""
gene_product_dict(model::EnzymeConstrainedModel, opt_model) =
    is_solved(opt_model) ?
    Dict(genes(model) .=> value.(opt_model[:x])[(length(model.columns)+1):end]) : nothing

"""
$(TYPEDSIGNATURES)

A pipe-able variant of [`gene_product_dict`](@ref).
"""
gene_product_dict(model::EnzymeConstrainedModel) = x -> gene_product_dict(model, x)

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

"""
$(TYPEDSIGNATURES)

Compute a "score" for picking the most viable isozyme for
[`make_simplified_enzyme_constrained_model`](@ref), based on maximum kcat divided by relative mass of
the isozyme. This is used because sMOMENT algorithm can not handle multiple
isozymes for one reaction.

# Note
This function does not take the direction of the reaction into account, i.e. the
maximum forward or reverse turnover number is used internally.
"""
simplified_enzyme_constrained_isozyme_speed(isozyme::Isozyme, gene_product_molar_mass) =
    max(isozyme.kcat_forward, isozyme.kcat_backward) /
    sum(count * gene_product_molar_mass(gene) for (gene, count) in isozyme.stoichiometry)

"""
$(TYPEDSIGNATURES)

A piping- and argmax-friendly overload of [`simplified_enzyme_constrained_isozyme_speed`](@ref).

# Example
```
gene_mass_function = gid -> 1.234

best_isozyme_for_simplified_enzyme_constrained = argmax(
    simplified_enzyme_constrained_isozyme_speed(gene_mass_function),
    my_isozyme_vector,
)
```
"""
simplified_enzyme_constrained_isozyme_speed(gene_product_molar_mass::Function) =
    isozyme -> simplified_enzyme_constrained_isozyme_speed(isozyme, gene_product_molar_mass)
