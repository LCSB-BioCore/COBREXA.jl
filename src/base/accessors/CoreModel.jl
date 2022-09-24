
"""
$(TYPEDSIGNATURES)

Get the reactions in a `CoreModel`.
"""
reactions(a::CoreModel)::Vector{String} = a.rxns

"""
$(TYPEDSIGNATURES)

Metabolites in a `CoreModel`.
"""
metabolites(a::CoreModel)::Vector{String} = a.mets

"""
$(TYPEDSIGNATURES)

`CoreModel` stoichiometry matrix.
"""
stoichiometry(a::CoreModel)::SparseMat = a.S

"""
$(TYPEDSIGNATURES)

`CoreModel` flux bounds.
"""
bounds(a::CoreModel)::Tuple{Vector{Float64},Vector{Float64}} = (a.xl, a.xu)

"""
$(TYPEDSIGNATURES)

`CoreModel` target flux balance.
"""
balance(a::CoreModel)::SparseVec = a.b

"""
$(TYPEDSIGNATURES)

`CoreModel` objective vector.
"""
objective(a::CoreModel)::SparseVec = a.c

"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::CoreModel, rid::String)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, first(indexin([rid], m.rxns))])...))

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction at index `ridx`.
"""
reaction_stoichiometry(m::CoreModel, ridx)::Dict{String,Float64} =
    Dict(m.mets[k] => v for (k, v) in zip(findnz(m.S[:, ridx])...))

"""
$(TYPEDSIGNATURES)

Retrieve a vector of all gene associations in a [`CoreModel`](@ref), in the
same order as `reactions(model)`.
"""
reaction_gene_association_vec(model::CoreModel)::Vector{Maybe{GeneAssociation}} = model.grrs

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction
index.
"""
reaction_gene_association(model::CoreModel, ridx::Int)::Maybe{GeneAssociation} =
    model.grrs[ridx]

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`CoreModel`](@ref) by reaction ID.
"""
reaction_gene_association(model::CoreModel, rid::String)::Maybe{GeneAssociation} =
    model.grrs[first(indexin([rid], model.rxns))]

"""
$(TYPEDSIGNATURES)

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
