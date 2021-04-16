
"""
    struct CoreCoupledModel <: MetabolicModel

The linear model with additional coupling constraints in the form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
mutable struct CoreCoupledModel <: MetabolicModel
    lm::CoreModel
    C::SparseMat
    cl::SparseVec
    cu::SparseVec

    function CoreCoupledModel(
        lm::MetabolicModel,
        C::M,
        cl::V,
        cu::V,
    ) where {V<:VecType,M<:MatType}

        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), n_reactions(lm)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))

        new(convert(CoreModel, lm), sparse(C), sparse(cl), sparse(cu))
    end
end

reactions(a::CoreCoupledModel) = reactions(a.lm)
metabolites(a::CoreCoupledModel) = metabolites(a.lm)
stoichiometry(a::CoreCoupledModel) = stoichiometry(a.lm)
bounds(a::CoreCoupledModel) = bounds(a.lm)
balance(a::CoreCoupledModel) = balance(a.lm)
objective(a::CoreCoupledModel) = objective(a.lm)

"""
    coupling(a::CoreCoupledModel)::SparseMat

Coupling constraint matrix for a `CoreCoupledModel`.
"""
function coupling(a::CoreCoupledModel)::SparseMat
    a.C
end

"""
    n_coupling_constraints(a::CoreCoupledModel)::Int

The number of coupling constraints in a `CoreCoupledModel`.
"""
function n_coupling_constraints(a::CoreCoupledModel)::Int
    return size(a.C, 1)
end

"""
    coupling_bounds(a::CoreCoupledModel)::Tuple{SparseVec,SparseVec}

Coupling bounds for a `CoreCoupledModel`.
"""
function coupling_bounds(a::CoreCoupledModel)::Tuple{SparseVec,SparseVec}
    (a.cl, a.cu)
end

"""
    Base.convert(::Type{CoreCoupledModel}, m::M) where {M <: MetabolicModel}

Make a `CoreCoupledModel` out of any compatible model type.
"""
function Base.convert(::Type{CoreCoupledModel}, m::M) where {M<:MetabolicModel}
    (cl, cu) = coupling_bounds(m)
    CoreCoupledModel(convert(CoreModel, m), coupling(m), cl, cu)
end
