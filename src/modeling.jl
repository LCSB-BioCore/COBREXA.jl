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
                      ub::AbstractFloat;
                      checkConsistency=false) where {V<:VT}
    return addReactions(m, reshape(s, (length(s), 1)), [c], [lb], [ub], checkConsistency=checkConsistency)
end

function addReactions(m::LinearModel,
                      s::V,
                      c::AbstractFloat,
                      lb::AbstractFloat,
                      ub::AbstractFloat,
                      names::String;
                      checkConsistency=false) where {V<:VT}
    return addReactions(m, reshape(s, (length(s), 1)), [c], [lb], [ub], [names], checkConsistency=checkConsistency)
end

function addReactions(m::LinearModel,
                      Sp::M,
                      c::V,
                      lb::V,
                      ub::V;
                      checkConsistency=false)  where {M<:MT,V<:VT}
    names = ["r$x" for x in length(m.rxns)+1:length(m.rxns)+length(ub)]
    return addReactions(m, Sp, c, lb, ub, names, checkConsistency=checkConsistency)
end

mutable struct ReactionStatus
           alreadyPresent::Bool
           index::Int
           info::String
end

function addReactions(m::LinearModel,
                      Sp::M,
                      c::V,
                      lb::V,
                      ub::V,
                      names::C;
                      checkConsistency=false) where {M<:MT,V<:VT,C<:ST}
    if checkConsistency
        newReactions = consistency(m, Sp, c, lb, ub, names)
    else
        newReactions = 1:length(names)
    end

    newS = hcat(m.S, Sp[:,newReactions])
    newc = vcat(m.c, c[newReactions])
    newlb = vcat(m.lb, lb[newReactions])
    newub = vcat(m.ub, ub[newReactions])
    newRxns = vcat(m.rxns, names[newReactions])
    return LinearModel(newS, m.b, newc, newlb, newub, newRxns, m.mets)
end


function consistency(m::LinearModel,
                          Sp::M,
                          c::V,
                          lb::V,
                          ub::V,
                          names::C) where {M<:MT,V<:VT,C<:ST}

    statuses = Array{ReactionStatus}(undef, length(names))
    for (i, name) in enumerate(names)
        index = findfirst(isequal(name), m.rxns)
        if isnothing(index)
            statuses[i] = ReactionStatus(false, 0, "new reaction")
        else
            statuses[i] = ReactionStatus(true, index, "reaction with the same name")
        end
    end

    newReactions = findall(x->x.alreadyPresent==true, statuses)
    return newReactions
end


"""
Returns the number of reactions in the LinearModel
"""
function nReactions(m::LinearModel)
    return length(m.rxns)
end

"""
Returns the number of metabolites in the LinearModel
"""
function nMetabolites(m::LinearModel)
    return length(m.mets)
end
