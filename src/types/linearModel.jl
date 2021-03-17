const VT = AbstractVector{Float64}
const Vt = SparseVector{Float64,Int64}
const MT = AbstractMatrix{Float64}
const Mt = SparseMatrixCSC{Float64,Int64}
const ST = AbstractVector{String}
const St = Vector{String}


"""
A linear optimization problem of the form:
```
min c^T x
s.t. S x = b
    cₗ ≤ C x ≤ cᵤ
    xₗ ≤ x ≤ xᵤ
```
"""
mutable struct LinearModel
    S::Mt
    b::Vt
    C::Mt
    cl::Vt
    cu::Vt
    c::Vt
    xl::Vt
    xu::Vt
    rxns::St
    mets::St

    function LinearModel(
        S::MX,
        b::VX,
        c::VX,
        xl::VX,
        xu::VX,
        rxns::SX,
        mets::SX,
    ) where {VX<:VT,MX<:MT,SX<:ST}

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
        S::MX,
        b::VX,
        C::MX2,
        cl::VX,
        cu::VX,
        c::VX,
        xl::VX,
        xu::VX,
        rxns::SX,
        mets::SX,
    ) where {VX<:VT,MX<:MT,MX2<:MT,SX<:ST}

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
