
"""
    struct CoreModel <: MetabolicModel

A "bare bones" core linear optimization problem of the form, with reaction and
metabolite names.
```
min c^T x
s.t. S x = b
      xₗ ≤ x ≤ xᵤ
```
"""
@with_repr mutable struct CoreModel <: MetabolicModel
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
    reactions(a::CoreModel)::Vector{String}

Get the reactions in a `CoreModel`.
"""
reactions(a::CoreModel)::Vector{String} = a.rxns

"""
    metabolites(a::CoreModel)::Vector{String}

Metabolites in a `CoreModel`.
"""
metabolites(a::CoreModel)::Vector{String} = a.mets

"""
    stoichiometry(a::CoreModel)::SparseMat

`CoreModel` stoichiometry matrix.
"""
stoichiometry(a::CoreModel)::SparseMat = a.S

"""
    bounds(a::CoreModel)::Tuple{Vector{Float64},Vector{Float64}}

`CoreModel` flux bounds.
"""
bounds(a::CoreModel)::Tuple{Vector{Float64},Vector{Float64}} = (a.xl, a.xu)

"""
    balance(a::CoreModel)::SparseVec

`CoreModel` target flux balance.
"""
balance(a::CoreModel)::SparseVec = a.b

"""
    objective(a::CoreModel)::SparseVec

`CoreModel` objective vector.
"""
objective(a::CoreModel)::SparseVec = a.c

"""
    genes(a::CoreModel)::Vector{String}

Collect all genes contained in the [`CoreModel`](@ref). The call is expensive
for large models, because the vector is not stored and instead gets rebuilt
each time this function is called.
"""
function genes(a::CoreModel)::Vector{String}
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
    reaction_stoichiometry(model::CoreModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::CoreModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

"""
    reaction_stoichiometry(model::CoreModel, ridx)::Dict{String, Float64}

Return the stoichiometry of reaction at index `ridx`.
"""
reaction_stoichiometry(m::CoreModel, ridx)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

"""
    reaction_gene_association_vec(model::CoreModel)::Vector{Maybe{GeneAssociation}}

Retrieve a vector of all gene associations in a [`CoreModel`](@ref), in the
same order as `reactions(model)`.
"""
reaction_gene_association_vec(model::CoreModel)::Vector{Maybe{GeneAssociation}} = model.grrs

"""
    reaction_gene_association(model::CoreModel, ridx::Int)::Maybe{GeneAssociation}

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction
index.
"""
reaction_gene_association(model::CoreModel, ridx::Int)::Maybe{GeneAssociation} =
    model.grrs[ridx]

"""
    reaction_gene_association(model::CoreModel, rid::String)::Maybe{GeneAssociation}

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction ID.
"""
reaction_gene_association(model::CoreModel, rid::String)::Maybe{GeneAssociation} =
    model.grrs[first(indexin([rid], model.rxns))]

"""
    Base.convert(::Type{CoreModel}, m::M) where {M <: MetabolicModel}

Make a `CoreModel` out of any compatible model type.
"""
function Base.convert(::Type{CoreModel}, m::M) where {M<:MetabolicModel}
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
