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
