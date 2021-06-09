
"""
    struct CoreModel <: MetabolicModel

A "bare bones" core linear optimization problem of the form, with reaction and
metabolite names.
```
min c^T x
s.t. S x = b
      xₗ ≤ x ≤ xᵤ
```
"""
mutable struct CoreModel <: MetabolicModel
    S::SparseMat
    b::SparseVec
    c::SparseVec
    xl::SparseVec
    xu::SparseVec
    rxns::Vector{String}
    mets::Vector{String}

    function CoreModel(
        S::M,
        b::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K,
    ) where {V<:VecType,M<:MatType,K<:StringVecType}

        all([length(b), length(mets)] .== size(S, 1)) ||
            throw(DimensionMismatch("inconsistent number of metabolites"))

        all([length(c), length(xl), length(xu), length(rxns)] .== size(S, 2)) ||
            throw(DimensionMismatch("inconsistent number of reactions"))

        new(sparse(S), sparse(b), sparse(c), sparse(xl), sparse(xu), rxns, mets)
    end
end

"""
    reactions(a::CoreModel)::Vector{String}

Get the reactions in a `CoreModel`.
"""
function reactions(a::CoreModel)::Vector{String}
    a.rxns
end

"""
    metabolites(a::CoreModel)::Vector{String}

Metabolites in a `CoreModel`.
"""
function metabolites(a::CoreModel)::Vector{String}
    a.mets
end

"""
    stoichiometry(a::CoreModel)::SparseMat

`CoreModel` stoichiometry matrix.
"""
function stoichiometry(a::CoreModel)::SparseMat
    a.S
end

"""
    bounds(a::CoreModel)::Tuple{SparseVec,SparseVec}

`CoreModel` flux bounds.
"""
function bounds(a::CoreModel)::Tuple{SparseVec,SparseVec}
    (a.xl, a.xu)
end

"""
    balance(a::CoreModel)::SparseVec

`CoreModel` target flux balance.
"""
function balance(a::CoreModel)::SparseVec
    a.b
end

"""
    objective(a::CoreModel)::SparseVec

`CoreModel` objective vector.
"""
function objective(a::CoreModel)::SparseVec
    a.c
end

"""
    reaction_stoichiometry(model::CoreModel, rxn_id::String)::Dict{String, Float64}

Return the reaction equation of reaction with id `rxn_id` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_stoichiometry(m::CoreModel, rxn_id::String)::Dict{String,Float64}
    Dict(
        m.mets[k] => v for
        (k, v) in zip(findnz(m.S[:, first(indexin([rxn_id], m.rxns))])...)
    )
end

"""
    reaction_stoichiometry(model::CoreModel, rxn_id::String)::Dict{String, Float64}

Return the reaction equation of reaction with id `rxn_ind` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_stoichiometry(m::CoreModel, rxn_ind::Int)::Dict{String,Float64}
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, rxn_ind])...))
end

"""
    Base.convert(::Type{CoreModel}, m::M) where {M <: MetabolicModel}

Make a `CoreModel` out of any compatible model type.
"""
function Base.convert(::Type{CoreModel}, m::M) where {M<:MetabolicModel}
    if typeof(m) == CoreModel
        return m
    end

    (xl, xu) = bounds(m)
    CoreModel(
        stoichiometry(m),
        balance(m),
        objective(m),
        xl,
        xu,
        reactions(m),
        metabolites(m),
    )
end
