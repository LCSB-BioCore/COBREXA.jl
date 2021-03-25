
"""
    abstract type MetabolicModel end

A helper supertype that wraps everything usable as a linear-like model for
COBREXA functions. If you want to use your own type, make it a subtype (so that
the functions typecheck) and add instances for the data accessor methods below.
"""
abstract type MetabolicModel end

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
const MatType = AbstractMatrix{Float64}
const VecType = AbstractVector{Float64}
const StringVecType = AbstractVector{String}

_missingImplError = (m, a) -> throw(MethodError(m, a))

"""
    reactions(a::LM)::Vector{String} where {LM<:MetabolicModel}

Return a vector of reaction identifiers in a model.
"""
function reactions(a::LM)::Vector{String} where {LM<:MetabolicModel}
    _missingImplError(reactions, (a,))
end

"""
    metabolites(a::LM)::Vector{String} where {LM<:MetabolicModel}

Return a vector of metabolite identifiers in a model.
"""
function metabolites(a::LM)::Vector{String} where {LM<:MetabolicModel}
    _missingImplError(metabolites, (a,))
end

"""
    nReactions(a::LM)::Int where {LM<:MetabolicModel}

Get the number of reactions in a model.
"""
function nReactions(a::LM)::Int where {LM<:MetabolicModel}
    length(reactions(a))
end

"""
    nMetabolites(a::LM)::Int where {LM<:MetabolicModel}

Get the number of metabolites in a model.
"""
function nMetabolites(a::LM)::Int where {LM<:MetabolicModel}
    length(metabolites(a))
end

"""
    stoichiometry(a::LM)::SparseMat where {LM<:MetabolicModel}

Get the sparse stoichiometry matrix of a model.
"""
function stoichiometry(a::LM)::SparseMat where {LM<:MetabolicModel}
    _missingImplError(stoichiometry, (a,))
end

"""
    bounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:MetabolicModel}

Get the lower and upper flux bounds of a model.
"""
function bounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:MetabolicModel}
    _missingImplError(bounds, (a,))
end

"""
    balance(a::LM)::SparseVec where {LM<:MetabolicModel}

Get the sparse balance vector of a model (ie. the `b` from `S x = b`).
"""
function balance(a::LM)::SparseVec where {LM<:MetabolicModel}
    _missingImplError(balance, (a,))
end

"""
    objective(a::LM)::SparseVec where {LM<:MetabolicModel}

Get the objective vector of a model.
"""
function objective(a::LM)::SparseVec where {LM<:MetabolicModel}
    _missingImplError(objective, (a,))
end

"""
    coupling(a::LM)::SparseMat where {LM<:MetabolicModel}

Get a matrix of coupling constraint definitions of a model.
"""
function coupling(a::LM)::SparseMat where {LM<:MetabolicModel}
    _missingImplError(coupling, (a,))
end

"""
    nCouplingConstraints(a::LM)::Int where {LM<:MetabolicModel}

Get the number of coupling constraints in a model.
"""
function nCouplingConstraints(a::LM)::Int where {LM<:MetabolicModel}
    size(coupling(a), 1)
end

"""
    couplingBounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:MetabolicModel}

Get the lower and upper bounds for each coupling bound in a model, as specified
by `coupling`.
"""
function couplingBounds(a::LM)::Tuple{SparseVec,SparseVec} where {LM<:MetabolicModel}
    _missingImplError(couplingBounds, (a,))
end
