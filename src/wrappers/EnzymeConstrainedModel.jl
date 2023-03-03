
"""
$(TYPEDEF)

A model with complex enzyme concentration and capacity bounds, as described in
*Sánchez, Benjamín J., et al. "Improving the phenotype predictions of a yeast
genome-scale metabolic model by incorporating enzymatic constraints." Molecular
systems biology 13.8 (2017): 935.*

Use [`make_enzyme_constrained_model`](@ref) or [`with_enzyme_constraints`](@ref)
to construct this kind of model.

The model wraps another "internal" model, and adds following modifications:
- enzymatic reactions with known enzyme information are split into multiple
  forward and reverse variants for each isozyme (affects the stoichiometric
  matrix),
- each enzyme reaction may have multiple variants per isozyme, thus the
  stoichiometric matrix will include all these virtual enzyme balances,
- reaction coupling is added to ensure the groups of isozyme reactions obey the
  global reaction flux bounds from the original model (affects the coupling),
- gene concentrations specified by each reaction and its gene product
  stoichiometry, can constrained by the user to reflect measurements, such as
  from mass spectrometry (affects the simple bounds),
- additional coupling is added to simulate total masses of different proteins
  grouped by type (e.g., membrane-bound and free-floating proteins), which can
  be again constrained by the user (this is slightly generalized from original
  GECKO algorithm, which only considers a single group of indiscernible
  proteins).

The structure contains fields `columns` that describe the contents of the
stoichiometry matrix columns, `coupling_row_reaction`,
`coupling_row_gene_product` and `coupling_row_mass_group` that describe
correspondence of the coupling rows to original model and determine the coupling
bounds (note: the coupling for gene product is actually added to stoichiometry,
not in [`coupling`](@ref)), and `inner`, which is the original wrapped model.
The `objective` of the model includes also the extra columns for individual
genes, as held by `coupling_row_gene_product`.

Implementation exposes the split reactions (available as `variables(model)`),
but retains the original "simple" reactions accessible by [`reactions`](@ref).
The related constraints are implemented using [`coupling`](@ref) and
[`coupling_bounds`](@ref).

To implement this wrapper for a model, the accessors
[`reaction_isozymes`](@ref), [`gene_product_lower_bound`](@ref),
[`gene_product_upper_bound](@ref), [`gene_product_molar_mass`](@ref), need to be
available. Additionally, the model needs to associate [`Isozyme`](@ref)s with
reactions. Reactions without enzymes, or those that should be ignored need to
return `nothing` when [`reaction_isozymes`](@ref) is called on them.

# Fields
$(TYPEDFIELDS)
"""
struct EnzymeConstrainedModel <: AbstractModelWrapper
    objective::SparseVec
    columns::Vector{EnzymeConstrainedReactionColumn}
    coupling_row_reaction::Vector{Int}
    coupling_row_gene_product::Vector{Tuple{Int,Tuple{Float64,Float64}}}
    coupling_row_mass_group::Vector{EnzymeConstrainedCapacity}

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::EnzymeConstrainedModel) = model.inner

function Accessors.stoichiometry(model::EnzymeConstrainedModel)
    irrevS = stoichiometry(model.inner) * enzyme_constrained_column_reactions(model)
    enzS = enzyme_constrained_gene_product_coupling(model)
    [
        irrevS spzeros(size(irrevS, 1), size(enzS, 1))
        -enzS I(size(enzS, 1))
    ]
end

Accessors.objective(model::EnzymeConstrainedModel) = model.objective

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

Accessors.n_variables(model::EnzymeConstrainedModel) =
    length(model.columns) + n_genes(model)

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
the original fluxes in the wrapped model.
"""
Accessors.reaction_variables(model::EnzymeConstrainedModel) =
    Accessors.Internal.make_mapping_dict(
        variables(model),
        reactions(model),
        reaction_variables_matrix(model),
    ) # TODO currently inefficient

"""
$(TYPEDSIGNATURES)

Get a mapping of enzyme concentration (on a mass basis, i.e. mass enzyme/mass
cell) variables to inner variables.
"""
Accessors.enzyme_variables(model::EnzymeConstrainedModel) =
    Dict(gid => Dict(gid => gene_product_molar_mass(model, gid)) for gid in genes(model)) # this is enough for all the semantics to work

"""
$(TYPEDSIGNATURES)

Get a mapping of enzyme groups to variables. See [`enzyme_variables`](@ref).
"""
function Accessors.enzyme_group_variables(model::EnzymeConstrainedModel)
    enz_ids = genes(model)
    Dict(
        grp.group_id => Dict(
            enz_ids[idx] => mm for
            (idx, mm) in zip(grp.gene_product_idxs, grp.gene_product_molar_masses)
        ) for grp in model.coupling_row_mass_group
    )
end

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

Accessors.n_coupling_constraints(model::EnzymeConstrainedModel) =
    n_coupling_constraints(model.inner) +
    length(model.coupling_row_reaction) +
    length(model.coupling_row_mass_group)

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

Accessors.balance(model::EnzymeConstrainedModel) =
    [balance(model.inner); spzeros(length(model.coupling_row_gene_product))]

Accessors.n_genes(model::EnzymeConstrainedModel) = length(model.coupling_row_gene_product)

Accessors.genes(model::EnzymeConstrainedModel) =
    genes(model.inner)[[idx for (idx, _) in model.coupling_row_gene_product]]

Accessors.metabolites(model::EnzymeConstrainedModel) =
    [metabolites(model.inner); genes(model) .* "#enzyme_constrained"]

Accessors.n_metabolites(model::EnzymeConstrainedModel) =
    n_metabolites(model.inner) + n_genes(model)
