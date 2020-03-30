"""
A bilevel linear optimization problem of the form:
```
min c^T * x
s.t. S x = b
     lb ≤ x ≤ ub
```
"""
mutable struct CobraLP
    S       ::AbstractMatrix
    b       ::Array{Float64,1}
    c       ::Array{Float64,1}
    lb      ::Array{Float64,1}
    ub      ::Array{Float64,1}
    rxns    ::Array{String,1}
    mets    ::Array{String,1}

    function CobraLP(
        S       ::AbstractMatrix,
        b       ::Array{Float64,1},
        c       ::Array{Float64,1},
        lb      ::Array{Float64,1},
        ub      ::Array{Float64,1},
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

function addReaction(m::CobraLP, s::Array{Float64,1}, b)
    newS = vcat(m.S, s')
    newb = vcat(m.b, b)
    return CobraLP(newS, newb, m.c, m.lb, m.ub, m.rxns, m.mets)
end

function addReactions(m::CobraLP, Sp::AbstractMatrix, b::Array{Float64,1})
    newS = vcat(m.S, Sp)
    newb = vcat(m.b, b)
    return CobraLP(newS, newb, m.c, m.lb, m.ub, m.rxns, m.mets)
end
