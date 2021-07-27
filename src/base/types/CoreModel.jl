
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
reactions(a::CoreModel)::Vector{String} = a.rxns

"""
    metabolites(a::CoreModel)::Vector{String}

Metabolites in a `CoreModel`.
"""
metabolites(a::CoreModel)::Vector{String} = a.mets

"""
    stoichiometry(a::CoreModel)::SparseMat

`CoreModel` stoichiometry matrix.
"""
stoichiometry(a::CoreModel)::SparseMat = a.S

"""
    bounds(a::CoreModel)::Tuple{SparseVec,SparseVec}

`CoreModel` flux bounds.
"""
bounds(a::CoreModel)::Tuple{SparseVec,SparseVec} = (a.xl, a.xu)

"""
    balance(a::CoreModel)::SparseVec

`CoreModel` target flux balance.
"""
balance(a::CoreModel)::SparseVec = a.b

"""
    objective(a::CoreModel)::SparseVec

`CoreModel` objective vector.
"""
objective(a::CoreModel)::SparseVec = a.c

"""
    reaction_stoichiometry(model::CoreModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::CoreModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

"""
    reaction_stoichiometry(model::CoreModel, ridx)::Dict{String, Float64}

Return the stoichiometry of reaction at index `ridx`.
"""
reaction_stoichiometry(m::CoreModel, ridx)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

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
