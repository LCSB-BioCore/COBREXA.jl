
"""
    struct CoupledLinearModel <: MetabolicModel

The linear model with additional coupling constraints in the form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
mutable struct CoupledLinearModel <: MetabolicModel
    lm::LinearModel
    C::SparseMat
    cl::SparseVec
    cu::SparseVec

    function CoupledLinearModel(
        lm::MetabolicModel,
        C::M,
        cl::V,
        cu::V,
    ) where {V<:VecType,M<:MatType}

        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), nReactions(lm)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))

        new(convert(LinearModel, lm), sparse(C), sparse(cl), sparse(cu))
    end
end

reactions(a::CoupledLinearModel) = reactions(a.lm)
metabolites(a::CoupledLinearModel) = metabolites(a.lm)
stoichiometry(a::CoupledLinearModel) = stoichiometry(a.lm)
bounds(a::CoupledLinearModel) = bounds(a.lm)
balance(a::CoupledLinearModel) = balance(a.lm)
objective(a::CoupledLinearModel) = objective(a.lm)

"""
    coupling(a::CoupledLinearModel)::SparseMat

Coupling constraint matrix for a `CoupledLinearModel`.
"""
function coupling(a::CoupledLinearModel)::SparseMat
    a.C
end

"""
    nCouplingConstraints(a::CoupledLinearModel)::Int

The number of coupling constraints in a `CoupledLinearModel`.
"""
function nCouplingConstraints(a::CoupledLinearModel)::Int
    return size(a.C, 1)
end

"""
    couplingBounds(a::CoupledLinearModel)::Tuple{SparseVec,SparseVec}

Coupling bounds for a `CoupledLinearModel`.
"""
function couplingBounds(a::CoupledLinearModel)::Tuple{SparseVec,SparseVec}
    (a.cl, a.cu)
end

"""
    Base.convert(::Type{CoupledLinearModel}, m::M) where {M <: MetabolicModel}

Make a `CoupledLinearModel` out of any compatible model type.
"""
function Base.convert(::Type{CoupledLinearModel}, m::M) where {M<:MetabolicModel}
    (cl, cu) = couplingBounds(m)
    CoupledLinearModel(convert(LinearModel, m), coupling(m), cl, cu)
end
