const VT = Union{Array{Float64,1}, SparseVector{Float64,Int64}}
const MT = Union{AbstractMatrix, SparseMatrixCSC{Float64,Int64}}
const ST = Union{Array{String,1}, SparseVector{String,Int64}}

"""
A linear optimization problem of the form:
```
min c^T x
s.t. S x = b
    cₗ ≤ C x ≤ cᵤ
    xₗ ≤ x ≤ xᵤ
```
"""
mutable struct LinearModel{M<:MT, V<:VT, K<:ST}
    S       ::M
    b       ::V
    C       ::M
    cl      ::V
    cu      ::V
    c       ::V
    xl      ::V
    xu      ::V
    rxns    ::K
    mets    ::K

    function LinearModel(
        S       ::M,
        b       ::V,
        c       ::V,
        xl      ::V,
        xu      ::V,
        rxns    ::K,
        mets    ::K) where {V<:VT,M<:MT,K<:ST}


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
        S       ::M1,
        b       ::V,
        C       ::M2,
        cl      ::V,
        cu      ::V,
        c       ::V,
        xl      ::V,
        xu      ::V,
        rxns    ::K,
        mets    ::K) where {V<:VT,M1<:MT,M2<:MT,K<:ST}

        checkInputDimensions(S, b, C, cl, cu, c, xl, xu, rxns, mets)

        sS = sparse(S)
        sb = sparse(b)
        sC = sparse(C)
        scl = sparse(cl)
        scu = sparse(cu)
        sc = sparse(c)
        sxl = sparse(xl)
        sxu = sparse(xu)

        new{SparseMatrixCSC{Float64,Int64}, SparseVector{Float64,Int64}, Array{String,1}}(sS, sb, sC, scl, scu, sc, sxl, sxu, rxns, mets)
    end
end
