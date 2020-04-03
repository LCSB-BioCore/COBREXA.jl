using SparseArrays

const VT = Union{Array{Float64,1}, SparseVector{Float64,Int64}}
const MT = Union{AbstractMatrix, SparseMatrixCSC{Float64,Int64}}
const ST = Union{Array{String,1}, SparseVector{String,Int64}}

"""
A linear optimization problem of the form:
```
min c^T * x
s.t. S x = b
     lb ≤ x ≤ ub
```
"""
mutable struct LinearModel{M<:MT, V<:VT, C<:ST}
    S       ::M
    b       ::V
    c       ::V
    lb      ::V
    ub      ::V
    rxns    ::C
    mets    ::C

    function LinearModel(
        S       ::M,
        b       ::V,
        c       ::V,
        lb      ::V,
        ub      ::V,
        rxns    ::C,
        mets    ::C) where {V<:VT,M<:MT,C<:ST}

        n_c = length(c)
        n_b = length(b)

        length(ub) == length(lb) || throw(DimensionMismatch("`lb` and `ub` don't have the same size"))
        n_c == length(lb) || throw(DimensionMismatch("`c` doesn't have the same size as `lb`"))

        size(S) == (n_b, n_c) || throw(DimensionMismatch("`S` shape doesn't match with `c` and `b`"))

        length(rxns) == n_c || throw(DimensionMismatch("`rxns` size doesn't match with `S`"))
        length(mets) == n_b || throw(DimensionMismatch("`mets` size doesn't match with `S`"))

        new{M,V,C}(S, b, c, lb, ub, rxns, mets)
    end
end


function addReactions(m::LinearModel,
                      s::V,
                      c::AbstractFloat,
                      lb::AbstractFloat,
                      ub::AbstractFloat) where {V<:VT}
    return addReactions(m, reshape(s, (length(s), 1)), [c], [lb], [ub])
end

function addReactions(m::LinearModel,
                      s::V,
                      c::AbstractFloat,
                      lb::AbstractFloat,
                      ub::AbstractFloat,
                      names::String) where {V<:VT}
    return addReactions(m, reshape(s, (length(s), 1)), [c], [lb], [ub], [names])
end

function addReactions(m::LinearModel,
                      Sp::M,
                      c::V,
                      lb::V,
                      ub::V)  where {M<:MT,V<:VT}
    names = ["r$x" for x in length(m.rxns)+1:length(m.rxns)+length(ub)]
    return addReactions(m, Sp, c, lb, ub, names)
end

function addReactions(m::LinearModel,
                      Sp::M,
                      c::V,
                      lb::V,
                      ub::V,
                      names::C) where {M<:MT,V<:VT,C<:ST}
    newS = hcat(m.S, Sp)
    newc = vcat(m.c, c)
    newlb = vcat(m.lb, lb)
    newub = vcat(m.ub, ub)
    newRxns = vcat(m.rxns, names)
    return LinearModel(newS, m.b, newc, newlb, newub, newRxns, m.mets)
end
