"""
    struct _gecko_column

A helper type for describing the contents of [`GeckoModel`](@ref)s.
"""
struct _gecko_column
    reaction_idx::Int
    isozyme_idx::Int
    direction::Int
    reaction_coupling_row::Int
    lb::Float64
    ub::Float64
    gene_product_coupling::Vector{Tuple{Int,Float64}}
end

"""
    struct GeckoModel <: ModelWrapper

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
coupling columns, `coupling_row_reaction`, `coupling_row_gene_product` and
`coupling_row_mass_group` that describe correspondence of the coupling rows to
original model and determine the coupling bounds, and `inner`, which is the
original wrapped model.

Implementation exposes the split reactions (available as `reactions(model)`),
but retains the original "simple" reactions accessible by [`fluxes`](@ref). All
constraints are implemented using [`coupling`](@ref) and
[`coupling_bounds`](@ref), i.e., all virtual metabolites described by GECKO are
purely virtual and do not occur in [`metabolites`](@ref).
"""
struct GeckoModel <: ModelWrapper
    objective::SparseVec
    columns::Vector{_gecko_column}
    coupling_row_reaction::Vector{Int}
    coupling_row_gene_product::Vector{Tuple{Int,Tuple{Float64,Float64}}}
    coupling_row_mass_group::Vector{Tuple{String, Vector{Int}, Vector{Float64}, Float64}}

    inner::MetabolicModel
end

unwrap_model(model::GeckoModel) = model.inner

"""
    stoichiometry(model::GeckoModel)

Return a stoichiometry of the [`GeckoModel`](@ref). The enzymatic reactions are
split into unidirectional forward and reverse ones, each of which may have
multiple variants per isozyme.
"""
function stoichiometry(model::GeckoModel)
    irrevS = stoichiometry(model.inner) * COBREXA._gecko_column_reactions(model)
    enzS = COBREXA._gecko_gene_product_coupling(model)
    [
        irrevS spzeros(size(irrevS, 1), size(enzS, 1))
        -enzS I(size(enzS, 1))
    ]
end

"""
    objective(model::GeckoModel)

Reconstruct an objective of the [`GeckoModel`](@ref).
"""
objective(model::GeckoModel) = model.objective

"""
    reactions(model::GeckoModel)

Returns the internal reactions in a [`GeckoModel`](@ref) (these may be split
to forward- and reverse-only parts with different isozyme indexes; reactions
IDs are mangled accordingly with suffixes), as well as the genes associated 
with enzymatic reactions.
"""
reactions(model::GeckoModel) = 
    let inner_reactions = reactions(model.inner)
        [
            _gecko_reaction_name(
                inner_reactions[col.reaction_idx],
                col.direction,
                col.isozyme_idx,
            ) for col in model.columns
        ]
    end

"""
    n_reactions(model::GeckoModel)

Returns the number of all irreversible reactions in `model` as well as the number of gene products 
that take part in enzymatic reactions.
"""
n_reactions(model::GeckoModel) = length(reactions(model))

"""
    bounds(model::GeckoModel)

Return variable bounds for [`GeckoModel`](@ref).
"""
function bounds(model::GeckoModel)
    lbs = [
        [col.lb for col in model.columns]
        [lb for (_, (lb, _)) in model.coupling_row_gene_product]
    ]
    ubs = [
        [col.ub for col in model.columns]
        [ub for (_, (_, ub)) in model.coupling_row_gene_product]
    ]
    (lbs, ubs)
end

"""
    reaction_flux(model::GeckoModel)

Get the mapping of the reaction rates in [`GeckoModel`](@ref) to the original
fluxes in the wrapped model.
"""
reaction_flux(model::GeckoModel) = _gecko_column_reactions(model)' * reaction_flux(model.inner) 

"""
    coupling(model::GeckoModel)

Return the coupling of [`GeckoModel`](@ref). That combines the coupling of the
wrapped model, coupling for split (arm) reactions, and the coupling for the total
enzyme capacity.
"""
function coupling(model::GeckoModel)
    innerC = coupling(model.inner) * _gecko_column_reactions(model)
    rxnC = _gecko_reaction_coupling(model)     
    enzcap = _gecko_mass_group_coupling(model)
    [
        innerC spzeros(size(innerC, 1), n_genes(model))
        rxnC spzeros(size(rxnC, 1), n_genes(model))
        spzeros(length(model.coupling_row_mass_group), length(model.columns)) enzcap
    ]    
end

"""
    n_coupling_constraints(model::GeckoModel)

Count the coupling constraints in [`GeckoModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::GeckoModel) =
    n_coupling_constraints(model.inner) +
    length(model.coupling_row_reaction) +
    length(model.coupling_row_mass_group)

"""
    coupling_bounds(model::GeckoModel)

The coupling bounds for [`GeckoModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
function coupling_bounds(model::GeckoModel)
    (iclb, icub) = coupling_bounds(model.inner)
    (ilb, iub) = bounds(model.inner)
    return (
        vcat(
            iclb,
            ilb[model.coupling_row_reaction],
            [0.0 for _ in model.coupling_row_mass_group],
        ),
        vcat(
            icub,
            iub[model.coupling_row_reaction],
            [ub for (_, _, _, ub) in model.coupling_row_mass_group],
        ),
    )
end

"""
    balance(model::GeckoModel)

Return the balance of the inner model, concatenated with a vector of 
zeros representing the enzyme balance of a [`GeckoModel`](@ref).
"""
balance(model::GeckoModel) = [balance(model.inner); spzeros(length(model.coupling_row_gene_product))]

"""
    n_genes(model::GeckoModel)

Return the number of genes that have enzymatic constraints associated with them.
"""
n_genes(model::GeckoModel) = length(genes(model))

"""
    genes(model::GeckoModel)

Return the gene ids of genes that have enzymatic constraints associated with them.  
"""
genes(model::GeckoModel) = genes(model.inner)[[idx for (idx, _) in model.coupling_row_gene_product]]

"""
    metabolites(model::GeckoModel)

Return the ids of all metabolites, both real and pseudo, for a [`GeckoModel`](@ref).
"""
metabolites(model::GeckoModel) = [metabolites(model.inner); genes(model).*"#supply"]

"""
    n_metabolites(model::GeckoModel)

Return the number of metabolites, both real and pseudo, for a [`GeckoModel`](@ref).
"""
n_metabolites(model::GeckoModel) = length(metabolites(model))