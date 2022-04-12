"""
    mutable struct SMomentModel <: MetabolicModel

Construct an enzyme capacity constrained model see `Bekiaris, Pavlos Stephanos,
and Steffen Klamt. "Automatic construction of metabolic models with enzyme
constraints." BMC bioinformatics, 2020.` for implementation details.

Note, `"§"` is reserved for internal use as a delimiter, no reaction id should
contain that character. Also note, SMOMENT assumes that each reaction only has a
single enzyme (one GRR) associated with it. It is required that a model be
modified to ensure that this condition is met. For ease-of-use,
[`remove_slow_isozymes!`](@ref) is supplied to effect this. Currently only
`modifications` that change attributes of the `optimizer` are supported.

# Fields
```
reaction_ids::Vector{String}
irrev_reaction_ids::Vector{String}
metabolites::Vector{String}
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
mutable struct SMomentModel <: MetabolicModel
    reaction_ids::Vector{String}
    irrev_reaction_ids::Vector{String}
    metabolites::Vector{String}
    c::SparseVec
    S::SparseMat
    b::SparseVec
    xl::SparseVec
    xu::SparseVec
end

"""
    stoichiometry(model::SMomentModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
stoichiometry(model::SMomentModel) = model.S

"""
    balance(model::SMomentModel)

Return stoichiometric balance.
"""
balance(model::SMomentModel) = model.b

"""
    objective(model::SMomentModel)

Return objective of `model`.
"""
objective(model::SMomentModel) = model.c

"""
    fluxes(model::SMomentModel)

Returns the reversible reactions in `model`. For 
the irreversible reactions, use [`reactions`][@ref].
"""
fluxes(model::SMomentModel) = model.reaction_ids

"""
    n_fluxes(model::SMomentModel)

Returns the number of reversible reactions in the model.
"""
n_fluxes(model::SMomentModel) = length(model.reaction_ids)

"""
    irreversible_reactions(model::SMomentModel)

Returns the irreversible reactions in `model`.
"""
reactions(model::SMomentModel) = model.irrev_reaction_ids

"""
    n_reactions(model::SMomentModel)

Returns the number of irreversible reactions in `model`.
"""
n_reactions(model::SMomentModel) = length(model.irrev_reaction_ids)

"""
    metabolites(model::SMomentModel)

Return the metabolites in `model`.
"""
metabolites(model::SMomentModel) = model.metabolites

"""
    n_metabolites(model::SMomentModel) = 

Return the number of metabolites in `model`.
"""
n_metabolites(model::SMomentModel) = length(metabolites(model))

"""
    bounds(model::SMomentModel)

Return variable bounds for `SMomentModel`.
"""
bounds(model::SMomentModel) = (model.xl, model.xu)

"""
    reaction_flux(model::MetabolicModel)

Helper function to get fluxes from optimization problem.
"""
function reaction_flux(model::SMomentModel)
    R = spzeros(n_fluxes(model), n_reactions(model) + 1)
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