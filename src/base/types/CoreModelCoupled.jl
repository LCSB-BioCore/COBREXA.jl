
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

"""
    reactions(a::CoreModelCoupled)

Extract reactions from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
reactions(a::CoreModelCoupled) = reactions(a.lm)

"""
    metabolites(a::CoreModelCoupled)

Extract metabolites from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
metabolites(a::CoreModelCoupled) = metabolites(a.lm)

"""
    stoichiometry(a::CoreModelCoupled)

Extract stoichiometry from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
stoichiometry(a::CoreModelCoupled) = stoichiometry(a.lm)

"""
    bounds(a::CoreModelCoupled)

Extract bounds from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
bounds(a::CoreModelCoupled) = bounds(a.lm)

"""
    balance(a::CoreModelCoupled)

Extract balance from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
balance(a::CoreModelCoupled) = balance(a.lm)

"""
    objective(a::CoreModelCoupled)

Extract objective from [`CoreModelCoupled`](@ref) (uses the internal
[`CoreModel`](@ref)).
"""
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
    reaction_equation(model::CoreModelCoupled, rxn_id::String)::Dict{String, Float64}

Return the reaction equation of reaction with id `rxn_id` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_equation(m::CoreModelCoupled, rxn_id::String)::Dict{String, Float64}
    reaction_equation(m.lm, rxn_id)
end

"""
    Base.convert(::Type{CoreModelCoupled}, mm::MetabolicModel)

Make a `CoreModelCoupled` out of any compatible model type.
"""
function Base.convert(::Type{CoreModelCoupled}, mm::MetabolicModel)
    if typeof(mm) == CoreModelCoupled
        return mm
    end

    (cl, cu) = coupling_bounds(mm)
    CoreModelCoupled(convert(CoreModel, mm), coupling(mm), cl, cu)
end
