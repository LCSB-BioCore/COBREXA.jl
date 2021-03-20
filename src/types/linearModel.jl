
"""
    abstract type AbstractCobraModel end

A helper supertype that wraps everything usable as a LinearModel for COBREXA
functions. If you want to use your own type, make it a subtype (so that the
functions typecheck) and add instances for the data accessor methods below.
"""
abstract type AbstractCobraModel end

const SpMtx = SparseMatrixCSC{Float64,Int}
const SpVec = SparseVector{Float64,Int}
const StrVec = Vector{String}
const MT = AbstractMatrix{Float64}
const VT = AbstractVector{Float64}
const ST = AbstractVector{String}

_missingImplError = (m, a) -> throw(MethodError(m, a))

function reactions(a::LM)::StrVec where {LM<:AbstractCobraModel}
    _missingImplError(reactions, (a,))
end

function metabolites(a::LM)::StrVec where {LM<:AbstractCobraModel}
    _missingImplError(metabolites, (a,))
end

function nReactions(a::LM)::Int where {LM<:AbstractCobraModel}
    length(reactions(a))
end

function nMetabolites(a::LM)::Int where {LM<:AbstractCobraModel}
    length(metabolites(a))
end

function stoichiometry(a::LM)::SpMtx where {LM<:AbstractCobraModel}
    _missingImplError(stoichiometry, (a,))
end

function bounds(a::LM)::Tuple{SpVec,SpVec} where {LM<:AbstractCobraModel}
    _missingImplError(bounds, (a,))
end

function balance(a::LM)::SpVec where {LM<:AbstractCobraModel}
    _missingImplError(balance, (a,))
end

function objective(a::LM)::SpVec where {LM<:AbstractCobraModel}
    _missingImplError(objective, (a,))
end

function coupling(a::LM)::SpMtx where {LM<:AbstractCobraModel}
    _missingImplError(coupling, (a,))
end

function nCouplingConstraints(a::LM)::Int where {LM<:AbstractCobraModel}
    size(coupling(a), 1)
end

function couplingBounds(a::LM)::Tuple{SpVec,SpVec} where {LM<:AbstractCobraModel}
    _missingImplError(couplingBounds, (a,))
end


"""
A concrete linear optimization problem of the form:
```
min c^T x
s.t. S x = b
    cₗ ≤ C x ≤ cᵤ
    xₗ ≤ x ≤ xᵤ
```

This is the default model supported by COBREXA model building functions.
"""
mutable struct LinearModel <: AbstractCobraModel
    S::SpMtx
    b::SpVec
    C::SpMtx
    cl::SpVec
    cu::SpVec
    c::SpVec
    xl::SpVec
    xu::SpVec
    rxns::StrVec
    mets::StrVec

    function LinearModel(
        S::M,
        b::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K,
    ) where {V<:VT,M<:MT,K<:ST}

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
    ) where {V<:VT,M1<:MT,M2<:MT,K<:ST}

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

function reactions(a::LinearModel)::StrVec
    a.rxns
end

function metabolites(a::LinearModel)::StrVec
    a.mets
end

function stoichiometry(a::LinearModel)::SpMtx
    a.S
end

function bounds(a::LinearModel)::Tuple{SpVec,SpVec}
    (a.xl, a.xu)
end

function balance(a::LinearModel)::SpVec
    a.b
end

function objective(a::LinearModel)::SpVec
    a.c
end

function coupling(a::LinearModel)::SpMtx
    a.C
end

function couplingBounds(a::LinearModel)::Tuple{SpVec,SpVec}
    (a.cl, a.cu)
end
