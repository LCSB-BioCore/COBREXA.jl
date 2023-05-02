
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

Accessors.variable_ids(a::MatrixModel)::Vector{String} = a.rxns

Accessors.Internal.@all_variables_are_reactions MatrixModel

Accessors.metabolite_ids(a::MatrixModel)::Vector{String} = a.mets

Accessors.stoichiometry(a::MatrixModel)::SparseMat = a.S

Accessors.variable_bounds(a::MatrixModel)::Tuple{Vector{Float64},Vector{Float64}} =
    (a.xl, a.xu)

Accessors.metabolite_bounds(a::MatrixModel)::SparseVec = a.b

Accessors.objective(a::MatrixModel)::SparseVec = a.c

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

Accessors.reaction_stoichiometry(m::MatrixModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

Accessors.reaction_stoichiometry(m::MatrixModel, ridx::Int)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

Accessors.reaction_gene_associations(
    model::MatrixModel,
    ridx::Int,
)::Maybe{GeneAssociationsDNF} = model.grrs[ridx]

Accessors.reaction_gene_associations(
    model::MatrixModel,
    rid::String,
)::Maybe{GeneAssociationsDNF} = model.grrs[first(indexin([rid], model.rxns))]

function Base.convert(::Type{MatrixModel}, m::M) where {M<:AbstractMetabolicModel}
    if typeof(m) == MatrixModel
        return m
    end

    (xl, xu) = variable_bounds(m)
    MatrixModel(
        stoichiometry(m),
        metabolite_bounds(m),
        objective(m),
        xl,
        xu,
        variable_ids(m),
        metabolite_ids(m),
        Vector{Maybe{GeneAssociationsDNF}}([
            reaction_gene_associations(m, id) for id in variable_ids(m)
        ]),
    )
end
