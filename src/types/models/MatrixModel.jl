
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
mutable struct MatrixModel <: AbstractMetabolicModel
    S::SparseMat
    b::SparseVec
    c::SparseVec
    xl::Vector{Float64}
    xu::Vector{Float64}
    rxns::Vector{String}
    mets::Vector{String}
    grrs::Vector{Maybe{GeneAssociationsDNF}}

    function MatrixModel(
        S::MatType,
        b::VecType,
        c::VecType,
        xl::VecType,
        xu::VecType,
        rxns::StringVecType,
        mets::StringVecType,
        grrs::Vector{Maybe{GeneAssociationsDNF}} = Vector{Maybe{GeneAssociationsDNF}}(
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

Get the reactions in a `MatrixModel`.
"""
Accessors.variables(a::MatrixModel)::Vector{String} = a.rxns

Accessors.Internal.@all_variables_are_reactions MatrixModel

"""
$(TYPEDSIGNATURES)

Metabolites in a `MatrixModel`.
"""
Accessors.metabolites(a::MatrixModel)::Vector{String} = a.mets

"""
$(TYPEDSIGNATURES)

`MatrixModel` stoichiometry matrix.
"""
Accessors.stoichiometry(a::MatrixModel)::SparseMat = a.S

"""
$(TYPEDSIGNATURES)

`MatrixModel` flux bounds.
"""
Accessors.bounds(a::MatrixModel)::Tuple{Vector{Float64},Vector{Float64}} = (a.xl, a.xu)

"""
$(TYPEDSIGNATURES)

`MatrixModel` target flux balance.
"""
Accessors.balance(a::MatrixModel)::SparseVec = a.b

"""
$(TYPEDSIGNATURES)

`MatrixModel` objective vector.
"""
Accessors.objective(a::MatrixModel)::SparseVec = a.c

"""
$(TYPEDSIGNATURES)

Collect all genes contained in the [`MatrixModel`](@ref). The call is expensive
for large models, because the vector is not stored and instead gets rebuilt
each time this function is called.
"""
function Accessors.genes(a::MatrixModel)::Vector{String}
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
Accessors.reaction_stoichiometry(m::MatrixModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction at index `ridx`.
"""
Accessors.reaction_stoichiometry(m::MatrixModel, ridx)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

"""
$(TYPEDSIGNATURES)

Retrieve the gene reaction associations from [`MatrixModel`](@ref) by reaction
index.
"""
Accessors.reaction_gene_associations(
    model::MatrixModel,
    ridx::Int,
)::Maybe{GeneAssociationsDNF} = model.grrs[ridx]

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`MatrixModel`](@ref) by reaction ID.
"""
Accessors.reaction_gene_associations(
    model::MatrixModel,
    rid::String,
)::Maybe{GeneAssociationsDNF} = model.grrs[first(indexin([rid], model.rxns))]

"""
$(TYPEDSIGNATURES)

Make a `MatrixModel` out of any compatible model type.
"""
function Base.convert(::Type{MatrixModel}, m::M) where {M<:AbstractMetabolicModel}
    if typeof(m) == MatrixModel
        return m
    end

    (xl, xu) = bounds(m)
    MatrixModel(
        stoichiometry(m),
        balance(m),
        objective(m),
        xl,
        xu,
        variables(m),
        metabolites(m),
        Vector{Maybe{GeneAssociationsDNF}}([
            reaction_gene_associations(m, id) for id in variables(m)
        ]),
    )
end
