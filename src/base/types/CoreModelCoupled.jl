
"""
$(TYPEDEF)

A matrix-based wrap that adds reaction coupling matrix to the inner model. A
flux `x` feasible in this model must satisfy:
```
    cₗ ≤ C x ≤ cᵤ
```

# Fields
$(TYPEDFIELDS)
"""
mutable struct CoreCoupling{M} <: ModelWrapper where {M<:MetabolicModel}
    lm::M
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}

    function CoreCoupling(
        lm::M,
        C::MatType,
        cl::VecType,
        cu::VecType,
    ) where {M<:MetabolicModel}
        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), size(lm.S, 2)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))
        # TODO if CoreCoupling is retained, the size call above needs to be
        # TODO replaced with n_reactions (which introduces ordering issues since
        # TODO accessors are defined after types.)
        new{M}(lm, sparse(C), collect(cl), collect(cu))
    end
end

"""
    const CoreModelCoupled = CoreCoupling{CoreModel}

A matrix-based linear model with additional coupling constraints in the form:
```
    cₗ ≤ C x ≤ cᵤ
```

Internally, the model is implemented using [`CoreCoupling`](@ref) that contains a single [`CoreModel`](@ref).
"""
const CoreModelCoupled = CoreCoupling{CoreModel}

CoreModelCoupled(lm::CoreModel, C::MatType, cl::VecType, cu::VecType) =
    CoreCoupling(lm, sparse(C), collect(cl), collect(cu))
