"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
