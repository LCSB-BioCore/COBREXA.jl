
"""
    struct LinearModel <: MetabolicModel

A concrete linear optimization problem of the form:
```
min c^T x
s.t. S x = b
      xₗ ≤ x ≤ xᵤ
```
"""
mutable struct LinearModel <: MetabolicModel
    S::SparseMat
    b::SparseVec
    c::SparseVec
    xl::SparseVec
    xu::SparseVec
    rxns::Vector{String}
    mets::Vector{String}

    function LinearModel(
        S::M,
        b::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K,
    ) where {V<:VecType,M<:MatType,K<:StringVecType}

        all([length(b), length(mets)] .== size(S, 1)) ||
            throw(DimensionMismatch("inconsistent number of reactions"))

        all([length(c), length(xl), length(xu), length(rxns)] .== size(S, 2)) ||
            throw(DimensionMismatch("inconsistent number of reactions"))

        new(sparse(S), sparse(b), sparse(c), sparse(xl), sparse(xu), rxns, mets)
    end
end

"""
    reactions(a::LinearModel)::Vector{String}

Get the reactions in a `LinearModel`.
"""
function reactions(a::LinearModel)::Vector{String}
    a.rxns
end

"""
    metabolites(a::LinearModel)::Vector{String}

Metabolites in a `LinearModel`.
"""
function metabolites(a::LinearModel)::Vector{String}
    a.mets
end

"""
    stoichiometry(a::LinearModel)::SparseMat

`LinearModel` stoichiometry matrix.
"""
function stoichiometry(a::LinearModel)::SparseMat
    a.S
end

"""
    bounds(a::LinearModel)::Tuple{SparseVec,SparseVec}

`LinearModel` flux bounds.
"""
function bounds(a::LinearModel)::Tuple{SparseVec,SparseVec}
    (a.xl, a.xu)
end

"""
    balance(a::LinearModel)::SparseVec

`LinearModel` target flux balance.
"""
function balance(a::LinearModel)::SparseVec
    a.b
end

"""
    objective(a::LinearModel)::SparseVec

`LinearModel` objective vector.
"""
function objective(a::LinearModel)::SparseVec
    a.c
end

"""
    coupling(a::LinearModel)::SparseMat

Coupling constraint matrix for a `LinearModel`, actually empty.
"""
function coupling(a::LinearModel)::SparseMat
    spzeros(0, nReactions(a))
end

"""
    couplingBounds(a::LinearModel)::Tuple{SparseVec,SparseVec}

Coupling bounds for a `LinearModel`, in fact always empty.
"""
function couplingBounds(a::LinearModel)::Tuple{SparseVec,SparseVec}
    (spzeros(0), spzeros(0))
end


"""
    Base.convert(::Type{LinearModel}, m::M) where {M <: MetabolicModel}

Make a `LinearModel` out of any compatible model type.
"""
function Base.convert(::Type{LinearModel}, m::M) where {M<:MetabolicModel}
    (xl, xu) = bounds(m)
    LinearModel(
        stoichiometry(m),
        balance(m),
        objective(m),
        xl,
        xu,
        reactions(m),
        metabolites(m),
    )
end
