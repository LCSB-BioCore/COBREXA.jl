
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
mutable struct MatrixCoupling{M} <: AbstractModelWrapper where {M<:AbstractMetabolicModel}
    lm::M
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}

    function MatrixCoupling(
        lm::M,
        C::MatType,
        cl::VecType,
        cu::VecType,
    ) where {M<:AbstractMetabolicModel}
        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), n_variables(lm)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))

        new{M}(lm, sparse(C), collect(cl), collect(cu))
    end
end

"""
$(TYPEDSIGNATURES)
"""
Accessors.unwrap_model(a::MatrixCoupling) = a.lm

"""
$(TYPEDSIGNATURES)
"""
Accessors.coupling(a::MatrixCoupling)::SparseMat = vcat(coupling(a.lm), a.C)

"""
$(TYPEDSIGNATURES)
"""
Accessors.n_coupling_constraints(a::MatrixCoupling)::Int =
    n_coupling_constraints(a.lm) + size(a.C, 1)

"""
$(TYPEDSIGNATURES)
"""
Accessors.coupling_bounds(a::MatrixCoupling)::Tuple{Vector{Float64},Vector{Float64}} =
    vcat.(coupling_bounds(a.lm), (a.cl, a.cu))

"""
$(TYPEDSIGNATURES)

Make a `MatrixCoupling` out of any compatible model type.
"""
function Base.convert(
    ::Type{MatrixCoupling{M}},
    mm::AbstractMetabolicModel;
    clone_coupling = true,
) where {M}
    if mm isa MatrixCoupling{M}
        mm
    elseif mm isa MatrixCoupling
        MatrixCoupling(convert(M, mm.lm), mm.C, mm.cl, mm.cu)
    elseif clone_coupling
        (cl, cu) = coupling_bounds(mm)
        MatrixCoupling(convert(M, mm), coupling(mm), cl, cu)
    else
        MatrixCoupling(convert(M, mm), spzeros(0, n_variables(mm)), spzeros(0), spzeros(0))
    end
end

"""
    const MatrixModelWithCoupling = MatrixCoupling{MatrixModel}

A matrix-based linear model with additional coupling constraints in the form:
```
    cₗ ≤ C x ≤ cᵤ
```

Internally, the model is implemented using [`MatrixCoupling`](@ref) that contains a single [`MatrixModel`](@ref).
"""
const MatrixModelWithCoupling = MatrixCoupling{MatrixModel}

MatrixModelWithCoupling(lm::MatrixModel, C::MatType, cl::VecType, cu::VecType) =
    MatrixCoupling(lm, sparse(C), collect(cl), collect(cu))

# these are special for MatrixModel-ish models
@inherit_model_methods MatrixModelWithCoupling (ridx::Int,) lm (ridx,) Accessors.reaction_stoichiometry
