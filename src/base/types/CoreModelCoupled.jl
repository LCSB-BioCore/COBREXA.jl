
"""
    struct CoreModelCoupled <: MetabolicModel

The linear model with additional coupling constraints in the form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
mutable struct CoreModelCoupled <: ModelWrapper
    lm::CoreModel
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}

    function CoreModelCoupled(lm::MetabolicModel, C::MatType, cl::VecType, cu::VecType)
        length(cu) == length(cl) ||
            throw(DimensionMismatch("`cl` and `cu` need to have the same size"))
        size(C) == (length(cu), n_reactions(lm)) ||
            throw(DimensionMismatch("wrong dimensions of `C`"))

        new(convert(CoreModel, lm), sparse(C), collect(cl), collect(cu))
    end
end

"""
    unwrap_model(a::CoreModelCoupled)

Get the internal [`CoreModel`](@ref) out of [`CoreModelCoupled`](@ref).
"""
unwrap_model(a::CoreModelCoupled) = a.lm

"""
    coupling(a::CoreModelCoupled)::SparseMat

Coupling constraint matrix for a `CoreModelCoupled`.
"""
coupling(a::CoreModelCoupled)::SparseMat = a.C

"""
    n_coupling_constraints(a::CoreModelCoupled)::Int

The number of coupling constraints in a `CoreModelCoupled`.
"""
n_coupling_constraints(a::CoreModelCoupled)::Int = size(a.C, 1)

"""
    coupling_bounds(a::CoreModelCoupled)::Tuple{Vector{Float64},Vector{Float64}}

Coupling bounds for a `CoreModelCoupled`.
"""
coupling_bounds(a::CoreModelCoupled)::Tuple{Vector{Float64},Vector{Float64}} = (a.cl, a.cu)

# these are special for CoreModel-ish models
@_inherit_model_methods CoreModelCoupled () lm () reaction_gene_association_vec
@_inherit_model_methods CoreModelCoupled (ridx::Int,) lm (ridx,) reaction_stoichiometry

"""
    Base.convert(::Type{CoreModelCoupled}, mm::MetabolicModel)

Make a `CoreModelCoupled` out of any compatible model type.
"""
function Base.convert(::Type{CoreModelCoupled}, mm::MetabolicModel)
    # TODO this might need a bit of rethinking and might be deprecated soon.
    # Eventually it seems to me that the coupling should be added as a
    # completely generic wrapper.
    if typeof(mm) == CoreModelCoupled
        return mm
    end

    (cl, cu) = coupling_bounds(mm)
    CoreModelCoupled(convert(CoreModel, mm), coupling(mm), cl, cu)
end
