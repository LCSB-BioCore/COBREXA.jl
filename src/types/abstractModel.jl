
"""
    abstract type AbstractCobraModel end

A helper supertype that wraps everything usable as a linear-like model for
COBREXA functions. If you want to use your own type, make it a subtype (so that
the functions typecheck) and add instances for the data accessor methods below.
"""
abstract type AbstractCobraModel end

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
const MatType = AbstractMatrix{Float64}
const VecType = AbstractVector{Float64}
const StringVecType = AbstractVector{String}

_missingImplError = (m, a) -> throw(MethodError(m, a))

function reactions(a::LM)::Vector{String} where {LM<:AbstractCobraModel}
    _missingImplError(reactions, (a,))
end

function metabolites(a::LM)::Vector{String} where {LM<:AbstractCobraModel}
    _missingImplError(metabolites, (a,))
end

function nReactions(a::LM)::Int where {LM<:AbstractCobraModel}
    length(reactions(a))
end

function nMetabolites(a::LM)::Int where {LM<:AbstractCobraModel}
    length(metabolites(a))
end

function stoichiometry(a::LM)::SparseMat where {LM<:AbstractCobraModel}
    _missingImplError(stoichiometry, (a,))
end

function bounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:AbstractCobraModel}
    _missingImplError(bounds, (a,))
end

function balance(a::LM)::SparseVec where {LM<:AbstractCobraModel}
    _missingImplError(balance, (a,))
end

function objective(a::LM)::SparseVec where {LM<:AbstractCobraModel}
    _missingImplError(objective, (a,))
end

function coupling(a::LM)::SparseMat where {LM<:AbstractCobraModel}
    _missingImplError(coupling, (a,))
end

function nCouplingConstraints(a::LM)::Int where {LM<:AbstractCobraModel}
    size(coupling(a), 1)
end

function couplingBounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:AbstractCobraModel}
    _missingImplError(couplingBounds, (a,))
end
