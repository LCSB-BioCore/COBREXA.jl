
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

_missing_impl_error = (m, a) -> throw(MethodError(m, a))

"""
    reactions(a::MetabolicModel)::Vector{String}

Return a vector of reaction identifiers in a model.
"""
function reactions(a::MetabolicModel)::Vector{String}
    _missing_impl_error(reactions, (a,))
end

"""
    metabolites(a::MetabolicModel)::Vector{String}

Return a vector of metabolite identifiers in a model.
"""
function metabolites(a::MetabolicModel)::Vector{String}
    _missing_impl_error(metabolites, (a,))
end

"""
    n_reactions(a::MetabolicModel)::Int

Get the number of reactions in a model.
"""
function n_reactions(a::MetabolicModel)::Int
    length(reactions(a))
end

"""
    n_metabolites(a::MetabolicModel)::Int

Get the number of metabolites in a model.
"""
function n_metabolites(a::MetabolicModel)::Int
    length(metabolites(a))
end

"""
    stoichiometry(a::MetabolicModel)::SparseMat

Get the sparse stoichiometry matrix of a model.
"""
function stoichiometry(a::MetabolicModel)::SparseMat
    _missing_impl_error(stoichiometry, (a,))
end

"""
    bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper flux bounds of a model.
"""
function bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}
    _missing_impl_error(bounds, (a,))
end

"""
    balance(a::MetabolicModel)::SparseVec

Get the sparse balance vector of a model (ie. the `b` from `S x = b`).
"""
function balance(a::MetabolicModel)::SparseVec
    _missing_impl_error(balance, (a,))
end

"""
    objective(a::MetabolicModel)::SparseVec

Get the objective vector of a model.
"""
function objective(a::MetabolicModel)::SparseVec
    _missing_impl_error(objective, (a,))
end

"""
    coupling(a::MetabolicModel)::SparseMat

Get a matrix of coupling constraint definitions of a model. By default, there
is no coupling in the models.
"""
function coupling(a::MetabolicModel)::SparseMat
    return spzeros(0, n_reactions(a))
end

"""
    n_coupling_constraints(a::MetabolicModel)::Int

Get the number of coupling constraints in a model.
"""
function n_coupling_constraints(a::MetabolicModel)::Int
    size(coupling(a), 1)
end

"""
    coupling_bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper bounds for each coupling bound in a model, as specified
by `coupling`. By default, the model does not have any coupling bounds.
"""
function coupling_bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}
    return (spzeros(0), spzeros(0))
end
