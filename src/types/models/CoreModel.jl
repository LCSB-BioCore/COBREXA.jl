
"""
$(TYPEDEF)

A "bare bones" core linear optimization problem of the form, with reaction and
metabolite names.
```
min c^T x
s.t. S x = b
      xₗ ≤ x ≤ xᵤ
```

# Fields
$(TYPEDFIELDS)
"""
mutable struct CoreModel <: AbstractMetabolicModel
    S::SparseMat
    b::SparseVec
    c::SparseVec
    xl::Vector{Float64}
    xu::Vector{Float64}
    rxns::Vector{String}
    mets::Vector{String}
    grrs::Vector{Maybe{GeneAssociation}}

    function CoreModel(
        S::MatType,
        b::VecType,
        c::VecType,
        xl::VecType,
        xu::VecType,
        rxns::StringVecType,
        mets::StringVecType,
        grrs::Vector{Maybe{GeneAssociation}} = Vector{Maybe{GeneAssociation}}(
            nothing,
            length(rxns),
        ),
    )
        all([length(b), length(mets)] .== size(S, 1)) ||
            throw(DimensionMismatch("inconsistent number of metabolites"))

        all(
            [length(c), length(xl), length(xu), length(rxns), length(grrs)] .== size(S, 2),
        ) || throw(DimensionMismatch("inconsistent number of reactions"))

        new(sparse(S), sparse(b), sparse(c), collect(xl), collect(xu), rxns, mets, grrs)
    end
end

"""
$(TYPEDSIGNATURES)

Get the reactions in a `CoreModel`.
"""
Accessors.reactions(a::CoreModel)::Vector{String} = a.rxns

"""
$(TYPEDSIGNATURES)

Metabolites in a `CoreModel`.
"""
Accessors.metabolites(a::CoreModel)::Vector{String} = a.mets

"""
$(TYPEDSIGNATURES)

`CoreModel` stoichiometry matrix.
"""
Accessors.stoichiometry(a::CoreModel)::SparseMat = a.S

"""
$(TYPEDSIGNATURES)

`CoreModel` flux bounds.
"""
Accessors.bounds(a::CoreModel)::Tuple{Vector{Float64},Vector{Float64}} = (a.xl, a.xu)

"""
$(TYPEDSIGNATURES)

`CoreModel` target flux balance.
"""
Accessors.balance(a::CoreModel)::SparseVec = a.b

"""
$(TYPEDSIGNATURES)

`CoreModel` objective vector.
"""
Accessors.objective(a::CoreModel)::SparseVec = a.c

"""
$(TYPEDSIGNATURES)

Collect all genes contained in the [`CoreModel`](@ref). The call is expensive
for large models, because the vector is not stored and instead gets rebuilt
each time this function is called.
"""
function Accessors.genes(a::CoreModel)::Vector{String}
    res = Set{String}()
    for grr in a.grrs
        isnothing(grr) && continue
        for gs in grr
            for g in gs
                push!(res, g)
            end
        end
    end
    sort(collect(res))
end

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
Accessors."""
Accessors.reaction_stoichiometry(m::CoreModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction at index `ridx`.
"""
Accessors.reaction_stoichiometry(m::CoreModel, ridx)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction
index.
"""
Accessors.reaction_gene_association(model::CoreModel, ridx::Int)::Maybe{GeneAssociation} =
    model.grrs[ridx]

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction ID.
"""
Accessors.reaction_gene_association(model::CoreModel, rid::String)::Maybe{GeneAssociation} =
    model.grrs[first(indexin([rid], model.rxns))]

"""
$(TYPEDSIGNATURES)

Make a `CoreModel` out of any compatible model type.
"""
function Base.convert(::Type{CoreModel}, m::M) where {M<:AbstractMetabolicModel}
    if typeof(m) == CoreModel
        return m
    end

    (xl, xu) = bounds(m)
    CoreModel(
        stoichiometry(m),
        balance(m),
        objective(m),
        xl,
        xu,
        reactions(m),
        metabolites(m),
        Vector{Maybe{GeneAssociation}}([
            reaction_gene_association(m, id) for id in reactions(m)
        ]),
    )
end
