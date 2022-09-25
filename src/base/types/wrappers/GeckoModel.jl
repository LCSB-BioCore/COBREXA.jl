"""
$(TYPEDEF)

A helper type for describing the contents of [`GeckoModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct _gecko_reaction_column
    reaction_idx::Int
    isozyme_idx::Int
    direction::Int
    reaction_coupling_row::Int
    lb::Float64
    ub::Float64
    gene_product_coupling::Vector{Tuple{Int,Float64}}
end

"""
$(TYPEDEF)

A helper struct that contains the gene product capacity terms organized by
the grouping type, e.g. metabolic or membrane groups etc.

# Fields
$(TYPEDFIELDS)
"""
struct _gecko_capacity
    group_id::String
    gene_product_idxs::Vector{Int}
    gene_product_molar_masses::Vector{Float64}
    group_upper_bound::Float64
end

"""
$(TYPEDEF)

A model with complex enzyme concentration and capacity bounds, as described in
*Sánchez, Benjamín J., et al. "Improving the phenotype predictions of a yeast
genome-scale metabolic model by incorporating enzymatic constraints." Molecular
systems biology 13.8 (2017): 935.*

Use [`make_gecko_model`](@ref) or [`with_gecko`](@ref) to construct this kind
of model.

The model wraps another "internal" model, and adds following modifications:
- enzymatic reactions with known enzyme information are split into multiple
  forward and reverse variants for each isozyme,
- reaction coupling is added to ensure the groups of isozyme reactions obey the
  global reaction flux bounds from the original model,
- gene concentrations specified by each reaction and its gene product stoichiometry,
  can constrained by the user to reflect measurements, such as
  from mass spectrometry,
- additional coupling is added to simulate total masses of different proteins
  grouped by type (e.g., membrane-bound and free-floating proteins), which can
  be again constrained by the user (this is slightly generalized from original
  GECKO algorithm, which only considers a single group of indiscernible
  proteins).

The structure contains fields `columns` that describe the contents of the
stoichiometry matrix columns, `coupling_row_reaction`,
`coupling_row_gene_product` and `coupling_row_mass_group` that describe
correspondence of the coupling rows to original model and determine the
coupling bounds (note: the coupling for gene product is actually added to
stoichiometry, not in [`coupling`](@ref)), and `inner`, which is the original
wrapped model. The `objective` of the model includes also the extra columns for
individual genes, as held by `coupling_row_gene_product`.

Implementation exposes the split reactions (available as `reactions(model)`),
but retains the original "simple" reactions accessible by [`fluxes`](@ref).
The related constraints are implemented using [`coupling`](@ref) and
[`coupling_bounds`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct GeckoModel <: ModelWrapper
    objective::SparseVec
    columns::Vector{_gecko_reaction_column}
    coupling_row_reaction::Vector{Int}
    coupling_row_gene_product::Vector{Tuple{Int,Tuple{Float64,Float64}}}
    coupling_row_mass_group::Vector{_gecko_capacity}

    inner::MetabolicModel
end


# Gecko

"""
$(TYPEDSIGNATURES)

Internal helper for systematically naming reactions in [`GeckoModel`](@ref).
"""
_gecko_reaction_name(original_name::String, direction::Int, isozyme_idx::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward#$isozyme_idx" :
    "$original_name#reverse#$isozyme_idx"

"""
$(TYPEDSIGNATURES)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
_gecko_reaction_column_reactions(model::GeckoModel) =
    _gecko_reaction_column_reactions(model.columns, model.inner)

"""
$(TYPEDSIGNATURES)

Helper method that doesn't require the whole [`GeckoModel`](@ref).
"""
_gecko_reaction_column_reactions(columns, inner) = sparse(
    [col.reaction_idx for col in columns],
    1:length(columns),
    [col.direction >= 0 ? 1 : -1 for col in columns],
    n_reactions(inner),
    length(columns),
)

"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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