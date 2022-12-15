
"""
$(TYPEDSIGNATURES)

Internal helper for systematically naming reactions in [`GeckoModel`](@ref).
"""
gecko_reaction_name(original_name::String, direction::Int, isozyme_idx::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward#$isozyme_idx" :
    "$original_name#reverse#$isozyme_idx"

"""
$(TYPEDSIGNATURES)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
gecko_column_reactions(model::GeckoModel) =
    gecko_column_reactions(model.columns, model.inner)

"""
$(TYPEDSIGNATURES)

Helper method that doesn't require the whole [`GeckoModel`](@ref).
"""
gecko_column_reactions(columns, inner) = sparse(
    [col.reaction_idx for col in columns],
    1:length(columns),
    [col.direction >= 0 ? 1 : -1 for col in columns],
    n_variables(inner),
    length(columns),
)

"""
$(TYPEDSIGNATURES)

Compute the part of the coupling for [`GeckoModel`](@ref) that limits the
"arm" reactions (which group the individual split unidirectional reactions).
"""
gecko_reaction_coupling(model::GeckoModel) =
    let tmp = [
            (col.reaction_coupling_row, i, col.direction) for
            (i, col) = enumerate(model.columns) if col.reaction_coupling_row != 0
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
$(TYPEDSIGNATURES)

Compute the part of the coupling for GeckoModel that limits the amount of each
kind of protein available.
"""
gecko_gene_product_coupling(model::GeckoModel) =
    let
        tmp = [
            (row, i, val) for (i, col) in enumerate(model.columns) for
            (row, val) in col.gene_product_coupling
        ]
        sparse(
            [row for (row, _, _) in tmp],
            [col for (_, col, _) in tmp],
            [val for (_, _, val) in tmp],
            length(model.coupling_row_gene_product),
            length(model.columns),
        )
    end

"""
$(TYPEDSIGNATURES)

Compute the part of the coupling for [`GeckoModel`](@ref) that limits the total
mass of each group of gene products.
"""
function gecko_mass_group_coupling(model::GeckoModel)
    tmp = [ # mm = molar mass, mg = mass group, i = row idx, j = col idx
        (i, j, mm) for (i, mg) in enumerate(model.coupling_row_mass_group) for
        (j, mm) in zip(mg.gene_product_idxs, mg.gene_product_molar_masses)
    ]
    sparse(
        [i for (i, _, _) in tmp],
        [j for (_, j, _) in tmp],
        [mm for (_, _, mm) in tmp],
        length(model.coupling_row_mass_group),
        n_genes(model),
    )
end
