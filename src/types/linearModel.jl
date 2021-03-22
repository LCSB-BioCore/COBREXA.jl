
"""
    struct LinearModel <: AbstractCobraModel

A concrete linear optimization problem of the form:
```
min c^T x
s.t. S x = b
    cₗ ≤ C x ≤ cᵤ
    xₗ ≤ x ≤ xᵤ
```
"""
mutable struct LinearModel <: AbstractCobraModel
    S::SparseMat
    b::SparseVec
    C::SparseMat
    cl::SparseVec
    cu::SparseVec
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

        sS = sparse(S)
        sb = sparse(b)
        sC = spzeros(0, length(c))
        scl = spzeros(0)
        scu = spzeros(0)
        sc = sparse(c)
        sxl = sparse(xl)
        sxu = sparse(xu)

        LinearModel(sS, sb, sC, scl, scu, sc, sxl, sxu, rxns, mets)
    end

    function LinearModel(
        S::M1,
        b::V,
        C::M2,
        cl::V,
        cu::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K,
    ) where {V<:VecType,M1<:MatType,M2<:MatType,K<:StringVecType}

        checkInputDimensions(S, b, C, cl, cu, c, xl, xu, rxns, mets)

        sS = sparse(S)
        sb = sparse(b)
        sC = sparse(C)
        scl = sparse(cl)
        scu = sparse(cu)
        sc = sparse(c)
        sxl = sparse(xl)
        sxu = sparse(xu)

        new(sS, sb, sC, scl, scu, sc, sxl, sxu, rxns, mets)
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

Coupling constraint matrix for a `LinearModel`.
"""
function coupling(a::LinearModel)::SparseMat
    a.C
end

"""
    couplingBounds(a::LinearModel)::Tuple{SparseVec,SparseVec}

Coupling bounds for a `LinearModel`.
"""
function couplingBounds(a::LinearModel)::Tuple{SparseVec,SparseVec}
    (a.cl, a.cu)
end
