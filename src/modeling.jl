using SparseArrays

"""
A linear optimization problem of the form:
```
min c^T * x
s.t. S x = b
     lb ≤ x ≤ ub
```
"""
mutable struct LinearModel
    S       ::Union{AbstractMatrix, SparseMatrixCSC{Float64,Int64}}
    b       ::Union{Array{Float64,1}, SparseVector{Float64,Int64}}
    c       ::Union{Array{Float64,1}, SparseVector{Float64,Int64}}
    lb      ::Union{Array{Float64,1}, SparseVector{Float64,Int64}}
    ub      ::Union{Array{Float64,1}, SparseVector{Float64,Int64}}
    rxns    ::Array{String,1}
    mets    ::Array{String,1}

    function LinearModel(
        S       ::Union{AbstractMatrix, SparseMatrixCSC{Float64,Int64}},
        b       ::Union{Array{Float64,1}, SparseVector{Float64,Int64}},
        c       ::Union{Array{Float64,1}, SparseVector{Float64,Int64}},
        lb      ::Union{Array{Float64,1}, SparseVector{Float64,Int64}},
        ub      ::Union{Array{Float64,1}, SparseVector{Float64,Int64}},
        rxns    ::Array{String,1},
        mets    ::Array{String,1})

        n_c = length(c)
        n_b = length(b)

        length(ub) == length(lb) || throw(DimensionMismatch("`lb` and `ub` don't have the same size"))
        n_c == length(lb) || throw(DimensionMismatch("`c` doesn't have the same size as `lb`"))

        size(S) == (n_b, n_c) || throw(DimensionMismatch("`S` shape doesn't match with `c` and `b`"))
        new(S, b, c, lb, ub, rxns, mets)
    end
end

function addReactions(m::LinearModel, s::Union{Array{Float64,1}, SparseVector{Float64,Int64}}, b)
    return addReactions(m, reshape(s, (1, length(s))), [b])
end

function addReactions(m::LinearModel, Sp::Union{AbstractMatrix, SparseMatrixCSC{Float64,Int64}}, b::Union{Array{Float64,1}, SparseVector{Float64,Int64}})
    newS = vcat(m.S, Sp)
    newb = vcat(m.b, b)
    return LinearModel(newS, newb, m.c, m.lb, m.ub, m.rxns, m.mets)
end
