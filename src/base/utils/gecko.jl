
"""
    _gecko_reaction_name(original_name::String, direction::Int)

Internal helper for systematically naming reactions in [`GeckoModel`](@ref).
"""
_gecko_reaction_name(original_name::String, direction::Int, isozyme_idx::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward#$isozyme_idx" :
    "$original_name#reverse#$isozyme_idx"

"""
    _gecko_reaction_column_reactions(model::GeckoModel)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
_gecko_reaction_column_reactions(model::GeckoModel) = sparse(
    [col.reaction_idx for col in model.columns],
    1:length(model.columns),
    [col.direction >= 0 ? 1 : -1 for col in model.columns],
    n_reactions(model.inner),
    length(model.columns),
)

"""
    _gecko_reaction_coupling(model::GeckoModel)

Compute the part of the coupling for [`GeckoModel`](@ref) that limits the
"arm" reactions (which group the individual split unidirectional reactions).
"""
_gecko_reaction_coupling(model::GeckoModel) =
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
    _gecko_gene_product_coupling(model::GeckoModel)

Compute the part of the coupling for GeckoModel that limits the amount of each
kind of protein available.
"""
_gecko_gene_product_coupling(model::GeckoModel) =
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
    _gecko_mass_group_coupling(model::GeckoModel)

Compute the part of the coupling for [`GeckoModel`](@ref) that limits the total
mass of each group of gene products.
"""
function _gecko_mass_group_coupling(model::GeckoModel)
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
