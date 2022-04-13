
"""
    _smoment_reaction_name(original_name::String, direction::Int)

Internal helper for systematically naming reactions in [`SMomentModel`](@ref).
"""
_smoment_reaction_name(original_name::String, direction::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward" : "$original_name#reverse"

"""
    _smoment_column_reactions(model::SMomentModel)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
_smoment_column_reactions(model::SMomentModel) = sparse(
    [col.reaction_idx for col in model.columns],
    1:length(model.columns),
    [col.direction >= 0 ? 1 : -1 for col in model.columns],
    n_reactions(model.inner),
    length(model.columns),
)

"""
    _smoment_reaction_coupling(model::SMomentModel)

Compute the part of the coupling for [`SMomentModel`](@ref) that limits the
"arm" reactions (which group the individual split unidirectional reactions).
"""
_smoment_reaction_coupling(model::SMomentModel) =
    let tmp = [
            (col.coupling_row, i, col.direction) for
            (i, col) = enumerate(model.columns) if col.direction != 0
        ]

        sparse(
            [row for (row, _, _) in tmp],
            [col for (_, col, _) in tmp],
            [val for (_, _, val) in tmp],
            length(model.coupling_row_reaction),
            length(model.columns),
        )
    end

"""
    smoment_isozyme_speed(isozyme::Isozyme, gene_product_molar_mass)

Compute a "score" for picking the most viable isozyme for
[`make_smoment_model`](@ref), based on maximum kcat divided by relative mass of
the isozyme. This is used because sMOMENT algorithm can not handle multiple
isozymes for one reaction.
"""
smoment_isozyme_speed(isozyme::Isozyme, gene_product_molar_mass) =
    max(isozyme.kcat_forward, isozyme.kcat_reverse) / sum(
        count * gene_product_molar_mass(gene) for
        (gene, count) in isozyme.gene_product_count
    )

"""
    smoment_isozyme_speed(gene_product_molar_mass::Function)

A piping- and argmax-friendly overload of [`smoment_isozyme_speed`](@ref).

# Example
```
gene_mass_function = gid -> 1.234

best_isozyme_for_smoment = argmax(
    smoment_isozyme_speed(gene_mass_function),
    my_isozyme_vector,
)
```
"""
smoment_isozyme_speed(gene_product_molar_mass::Function) =
    isozyme -> smoment_isozyme_speed(isozyme, gene_product_molar_mass)
