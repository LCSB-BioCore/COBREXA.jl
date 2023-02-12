"""
$(TYPEDEF)

A helper type for describing the contents of [`EnzymeConstrainedModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct _EnzymeConstrainedReactionColumn
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
struct _EnzymeConstrainedCapacity
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

Use [`make_enzyme_constrained_model`](@ref) or [`with_enzyme_constraints`](@ref) to construct this kind
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

Implementation exposes the split reactions (available as `variables(model)`),
but retains the original "simple" reactions accessible by [`reactions`](@ref).
The related constraints are implemented using [`coupling`](@ref) and
[`coupling_bounds`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct EnzymeConstrainedModel <: AbstractModelWrapper
    objective::SparseVec
    columns::Vector{_EnzymeConstrainedReactionColumn}
    coupling_row_reaction::Vector{Int}
    coupling_row_gene_product::Vector{Tuple{Int,Tuple{Float64,Float64}}}
    coupling_row_mass_group::Vector{_EnzymeConstrainedCapacity}

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::EnzymeConstrainedModel) = model.inner

"""
$(TYPEDSIGNATURES)

Return a stoichiometry of the [`EnzymeConstrainedModel`](@ref). The enzymatic reactions are
split into unidirectional forward and reverse ones, each of which may have
multiple variants per isozyme.
"""
function Accessors.stoichiometry(model::EnzymeConstrainedModel)
    irrevS = stoichiometry(model.inner) * enzyme_constrained_column_reactions(model)
    enzS = enzyme_constrained_gene_product_coupling(model)
    [
        irrevS spzeros(size(irrevS, 1), size(enzS, 1))
        -enzS I(size(enzS, 1))
    ]
end

"""
$(TYPEDSIGNATURES)

Return the objective of the [`EnzymeConstrainedModel`](@ref). Note, the objective is with
respect to the internal variables, i.e. [`variables(model)`](@ref), which are
the unidirectional reactions and the genes involved in enzymatic reactions that
have kinetic data.
"""
Accessors.objective(model::EnzymeConstrainedModel) = model.objective

"""
$(TYPEDSIGNATURES)

Returns the internal reactions in a [`EnzymeConstrainedModel`](@ref) (these may be split
to forward- and reverse-only parts with different isozyme indexes; reactions
IDs are mangled accordingly with suffixes).
"""
function Accessors.variables(model::EnzymeConstrainedModel)
    inner_reactions = variables(model.inner)
    mangled_reactions = [
        enzyme_constrained_reaction_name(
            inner_reactions[col.reaction_idx],
            col.direction,
            col.isozyme_idx,
        ) for col in model.columns
    ]
    [mangled_reactions; genes(model)]
end

"""
$(TYPEDSIGNATURES)

Returns the number of all irreversible reactions in `model` as well as the
number of gene products that take part in enzymatic reactions.
"""
Accessors.n_variables(model::EnzymeConstrainedModel) =
    length(model.columns) + n_genes(model)

"""
$(TYPEDSIGNATURES)

Return variable bounds for [`EnzymeConstrainedModel`](@ref).
"""
function Accessors.bounds(model::EnzymeConstrainedModel)
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
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`EnzymeConstrainedModel`](@ref) to
the original fluxes in the wrapped model (as a matrix).
"""
function Accessors.reaction_variables_matrix(model::EnzymeConstrainedModel)
    rxnmat =
        enzyme_constrained_column_reactions(model)' * reaction_variables_matrix(model.inner)
    [
        rxnmat
        spzeros(n_genes(model), size(rxnmat, 2))
    ]
end

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`EnzymeConstrainedModel`](@ref) to
the original fluxes in the wrapped model
"""
Accessors.reaction_variables(model::EnzymeConstrainedModel) =
    Accessors.Internal.make_mapping_dict(
        variables(model),
        reactions(model),
        reaction_variables_matrix(model),
    ) # TODO currently inefficient

"""
$(TYPEDSIGNATURES)

Get a mapping of enzyme variables to variables -- for enzyme constrained models,
this is just a direct mapping.
"""
Accessors.enzyme_variables(model::EnzymeConstrainedModel) = Dict(
    gid => Dict(gid => 1.0) for gid in genes(model)
) # this is enough for all the semantics to work

"""
$(TYPEDSIGNATURES)

Return the coupling of [`EnzymeConstrainedModel`](@ref). That combines the coupling of the
wrapped model, coupling for split (arm) reactions, and the coupling for the total
enzyme capacity.
"""
function Accessors.coupling(model::EnzymeConstrainedModel)
    innerC = coupling(model.inner) * enzyme_constrained_column_reactions(model)
    rxnC = enzyme_constrained_reaction_coupling(model)
    enzcap = enzyme_constrained_mass_group_coupling(model)
    [
        innerC spzeros(size(innerC, 1), n_genes(model))
        rxnC spzeros(size(rxnC, 1), n_genes(model))
        spzeros(length(model.coupling_row_mass_group), length(model.columns)) enzcap
    ]
end

"""
$(TYPEDSIGNATURES)

Count the coupling constraints in [`EnzymeConstrainedModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
Accessors.n_coupling_constraints(model::EnzymeConstrainedModel) =
    n_coupling_constraints(model.inner) +
    length(model.coupling_row_reaction) +
    length(model.coupling_row_mass_group)

"""
$(TYPEDSIGNATURES)

The coupling bounds for [`EnzymeConstrainedModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
function Accessors.coupling_bounds(model::EnzymeConstrainedModel)
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
            [grp.group_upper_bound for grp in model.coupling_row_mass_group],
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Return the balance of the reactions in the inner model, concatenated with a vector of
zeros representing the enzyme balance of a [`EnzymeConstrainedModel`](@ref).
"""
Accessors.balance(model::EnzymeConstrainedModel) =
    [balance(model.inner); spzeros(length(model.coupling_row_gene_product))]

"""
$(TYPEDSIGNATURES)

Return the number of genes that have enzymatic constraints associated with them.
"""
Accessors.n_genes(model::EnzymeConstrainedModel) = length(model.coupling_row_gene_product)

"""
$(TYPEDSIGNATURES)

Return the gene ids of genes that have enzymatic constraints associated with them.
"""
Accessors.genes(model::EnzymeConstrainedModel) =
    genes(model.inner)[[idx for (idx, _) in model.coupling_row_gene_product]]

"""
$(TYPEDSIGNATURES)

Return the ids of all metabolites, both real and pseudo, for a [`EnzymeConstrainedModel`](@ref).
"""
Accessors.metabolites(model::EnzymeConstrainedModel) =
    [metabolites(model.inner); genes(model) .* "#enzyme_constrained"]

"""
$(TYPEDSIGNATURES)

Return the number of metabolites, both real and pseudo, for a [`EnzymeConstrainedModel`](@ref).
"""
Accessors.n_metabolites(model::EnzymeConstrainedModel) =
    n_metabolites(model.inner) + n_genes(model)
