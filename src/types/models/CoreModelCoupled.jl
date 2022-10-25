
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
mutable struct CoreCoupling{M} <: ModelWrapper where {M<:AbstractMetabolicModel}
    lm::M
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}

    function CoreCoupling(
        lm::M,
        C::MatType,
        cl::VecType,
        cu::VecType,
    ) where {M<:AbstractMetabolicModel}
        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), n_reactions(lm)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))

        new{M}(lm, sparse(C), collect(cl), collect(cu))
    end
end

"""
$(TYPEDSIGNATURES)

Get the internal [`CoreModel`](@ref) out of [`CoreCoupling`](@ref).
"""
Accessors.unwrap_model(a::CoreCoupling) = a.lm

"""
$(TYPEDSIGNATURES)

Coupling constraint matrix for a `CoreCoupling`.
"""
Accessors.coupling(a::CoreCoupling)::SparseMat = vcat(coupling(a.lm), a.C)

"""
$(TYPEDSIGNATURES)

The number of coupling constraints in a `CoreCoupling`.
"""
Accessors.n_coupling_constraints(a::CoreCoupling)::Int =
    n_coupling_constraints(a.lm) + size(a.C, 1)

"""
$(TYPEDSIGNATURES)

Coupling bounds for a `CoreCoupling`.
"""
Accessors.coupling_bounds(a::CoreCoupling)::Tuple{Vector{Float64},Vector{Float64}} =
    vcat.(coupling_bounds(a.lm), (a.cl, a.cu))

"""
$(TYPEDSIGNATURES)

Make a `CoreCoupling` out of any compatible model type.
"""
function Base.convert(
    ::Type{CoreCoupling{M}},
    mm::AbstractMetabolicModel;
    clone_coupling = true,
) where {M}
    if mm isa CoreCoupling{M}
        mm
    elseif mm isa CoreCoupling
        CoreCoupling(convert(M, mm.lm), mm.C, mm.cl, mm.cu)
    elseif clone_coupling
        (cl, cu) = coupling_bounds(mm)
        CoreCoupling(convert(M, mm), coupling(mm), cl, cu)
    else
        CoreCoupling(convert(M, mm), spzeros(0, n_reactions(mm)), spzeros(0), spzeros(0))
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

# these are special for CoreModel-ish models
@inherit_model_methods CoreModelCoupled (ridx::Int,) lm (ridx,) Accessors.reaction_stoichiometry
