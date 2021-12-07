
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

"""
    reaction_stoichiometry(model::CoreModelCoupled, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::CoreModelCoupled, rid::String) = reaction_stoichiometry(m.lm, rid)

"""
    reaction_stoichiometry(model::CoreModelCoupled, ridx)::Dict{String, Float64}

Return the stoichiometry of reaction at index `ridx`.
"""
function reaction_stoichiometry(m::CoreModelCoupled, ridx)::Dict{String,Float64}
    reaction_stoichiometry(m.lm, ridx)
end

"""
    reaction_gene_association_vec(model::CoreModelCoupled)::Vector{Maybe{GeneAssociation}}

Retrieve a vector of gene associations in a [`CoreModelCoupled`](@ref), in the
same order as `reactions(model)`.
"""
reaction_gene_association_vec(model::CoreModelCoupled)::Vector{Maybe{GeneAssociation}} =
    reaction_gene_association_vec(model.lm)

"""
    reaction_gene_association(model::CoreModelCoupled, ridx::Int)::Maybe{GeneAssociation}

Retrieve the [`GeneAssociation`](@ref) from [`CoreModelCoupled`](@ref) by reaction
index.
"""
reaction_gene_association(model::CoreModelCoupled, ridx::Int)::Maybe{GeneAssociation} =
    reaction_gene_association(model.lm, ridx)

"""
    reaction_gene_association(model::CoreModelCoupled, rid::String)::Maybe{GeneAssociation}

Retrieve the [`GeneAssociation`](@ref) from [`CoreModelCoupled`](@ref) by reaction ID.
"""
reaction_gene_association(model::CoreModelCoupled, rid::String)::Maybe{GeneAssociation} =
    reaction_gene_association(model.lm, rid)

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
