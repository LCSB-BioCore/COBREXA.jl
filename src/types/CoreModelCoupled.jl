
"""
    struct CoreModelCoupled <: MetabolicModel

The linear model with additional coupling constraints in the form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
mutable struct CoreModelCoupled <: MetabolicModel
    lm::CoreModel
    C::SparseMat
    cl::SparseVec
    cu::SparseVec

    function CoreModelCoupled(
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

reactions(a::CoreModelCoupled) = reactions(a.lm)
metabolites(a::CoreModelCoupled) = metabolites(a.lm)
stoichiometry(a::CoreModelCoupled) = stoichiometry(a.lm)
bounds(a::CoreModelCoupled) = bounds(a.lm)
balance(a::CoreModelCoupled) = balance(a.lm)
objective(a::CoreModelCoupled) = objective(a.lm)

"""
    coupling(a::CoreModelCoupled)::SparseMat

Coupling constraint matrix for a `CoreModelCoupled`.
"""
function coupling(a::CoreModelCoupled)::SparseMat
    a.C
end

"""
    n_coupling_constraints(a::CoreModelCoupled)::Int

The number of coupling constraints in a `CoreModelCoupled`.
"""
function n_coupling_constraints(a::CoreModelCoupled)::Int
    return size(a.C, 1)
end

"""
    coupling_bounds(a::CoreModelCoupled)::Tuple{SparseVec,SparseVec}

Coupling bounds for a `CoreModelCoupled`.
"""
function coupling_bounds(a::CoreModelCoupled)::Tuple{SparseVec,SparseVec}
    (a.cl, a.cu)
end

"""
    Base.convert(::Type{CoreModelCoupled}, m::M) where {M <: MetabolicModel}

Make a `CoreModelCoupled` out of any compatible model type.
"""
function Base.convert(::Type{CoreModelCoupled}, m::M) where {M<:MetabolicModel}
    (cl, cu) = coupling_bounds(m)
    CoreModelCoupled(convert(CoreModel, m), coupling(m), cl, cu)
end
