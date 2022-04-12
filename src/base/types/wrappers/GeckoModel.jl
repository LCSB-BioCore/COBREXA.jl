"""
    mutable struct GeckoModel <: MetabolicModel

A model that incorporates enzyme capacity and kinetic constraints via the GECKO
formulation. See `Sánchez, Benjamín J., et al. "Improving the phenotype
predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
constraints." Molecular systems biology, 2017.` for implementation details.

Note, since the model uses irreversible reactions internally, `"§FOR"` (for the
forward direction) and `"§REV"` (for the reverse direction) is appended to each
reaction internally. Hence, `"§"` is reserved for internal use as a delimiter,
no reaction id should contain this character. 

To actually run GECKO, call [`flux_balance_analysis`](@ref) on a `GeckoModel`.

# Fields
```
reaction_ids::Vector{String}
irrev_reaction_ids::Vector{String}
metabolites::Vector{String}
gene_ids::Vector{String}
c::SparseVec
S::SparseMat
b::SparseVec
xl::SparseVec
xu::SparseVec
C::SparseMat
cl::Vector{Float64}
cu::Vector{Float64}
```
"""
mutable struct GeckoModel <: MetabolicModel
    reaction_ids::Vector{String}
    irrev_reaction_ids::Vector{String}
    metabolites::Vector{String}
    gene_ids::Vector{String}

    # gecko  
    c::SparseVec
    S::SparseMat
    b::SparseVec
    xl::SparseVec
    xu::SparseVec

    # enzyme capacity constraints
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}
end

"""
    stoichiometry(model::GeckoModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
stoichiometry(model::GeckoModel) = model.S

"""
    balance(model::GeckoModel)

Return stoichiometric balance.
"""
balance(model::GeckoModel) = model.b

"""
    objective(model::GeckoModel)

Return objective of `model`.
"""
objective(model::GeckoModel) = model.c

"""
    fluxes(model::GeckoModel)

Returns the reversible reactions in `model`. For 
the irreversible reactions, use [`reactions`][@ref].
"""
fluxes(model::GeckoModel) = model.reaction_ids

"""
    n_reactions(model::GeckoModel)

Returns the number of reversible reactions in the model.
"""
n_fluxes(model::GeckoModel) = length(model.reaction_ids)

"""
    reactions(model::GeckoModel)

Returns the irreversible reactions in `model`.
"""
reactions(model::GeckoModel) = model.irrev_reaction_ids

"""
    reactions(model::GeckoModel)

Returns the number of all irreversible reactions in `model`.
"""
n_reactions(model::GeckoModel) = length(model.irrev_reaction_ids)

"""
    genes(model::GeckoModel)

Returns the genes (proteins) in the order as they appear as variables in the
model.
"""
genes(model::GeckoModel) = model.gene_ids

"""
    n_genes(model::GeckoModel)

Returns the number of genes in the model.
"""
n_genes(model::GeckoModel) = length(model.gene_ids)

"""
    metabolites(model::GeckoModel)

Return the metabolites in `model`.
"""
metabolites(model::GeckoModel) = model.metabolites

"""
    n_metabolites(model::GeckoModel) = 

Return the number of metabolites in `model`.
"""
n_metabolites(model::GeckoModel) = length(metabolites(model))

"""
    bounds(model::GeckoModel)

Return variable bounds for `GeckoModel`.
"""
bounds(model::GeckoModel) = (model.xl, model.xu)

"""
    coupling(model::GeckoModel)

Coupling constraint matrix for a `GeckoModel`.
"""
coupling(model::GeckoModel) = model.C

"""
    coupling_bounds(model::GeckoModel)

Coupling bounds for a `GeckoModel`.
"""
coupling_bounds(model::GeckoModel) = (model.cl, model.cu)

"""
    reaction_flux(model::MetabolicModel)

Helper function to get fluxes from optimization problem.
"""
function reaction_flux(model::GeckoModel)
    R = spzeros(n_fluxes(model), n_genes(model) + n_reactions(model))
    for (i, rid) in enumerate(fluxes(model))
        for_idx = findfirst(
            x -> x == rid * "§ARM§FOR" || x == rid * "§FOR",
            model.irrev_reaction_ids,
        )
        rev_idx = findfirst(
            x -> x == rid * "§ARM§REV" || x == rid * "§REV",
            model.irrev_reaction_ids,
        )
        !isnothing(for_idx) && (R[i, for_idx] = 1.0)
        !isnothing(rev_idx) && (R[i, rev_idx] = -1.0)
    end
    return R'
end
