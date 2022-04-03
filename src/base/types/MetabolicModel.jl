
#
# IMPORTANT
#
# This file provides a list of "officially supported" accessors that should
# work with all subtypes of [`MetabolicModel`](@ref). Keep this synced with the
# automatically derived methods for [`ModelWrapper`](@ref).
#

_missing_impl_error(m, a) = throw(MethodError(m, a))

"""
    reactions(a::MetabolicModel)::Vector{String}

Return a vector of reaction identifiers in a model.
"""
function reactions(a::MetabolicModel)::Vector{String}
    _missing_impl_error(reactions, (a,))
end

"""
    metabolites(a::MetabolicModel)::Vector{String}

Return a vector of metabolite identifiers in a model.
"""
function metabolites(a::MetabolicModel)::Vector{String}
    _missing_impl_error(metabolites, (a,))
end

"""
    n_reactions(a::MetabolicModel)::Int

Get the number of reactions in a model.
"""
function n_reactions(a::MetabolicModel)::Int
    length(reactions(a))
end

"""
    n_metabolites(a::MetabolicModel)::Int

Get the number of metabolites in a model.
"""
function n_metabolites(a::MetabolicModel)::Int
    length(metabolites(a))
end

"""
    stoichiometry(a::MetabolicModel)::SparseMat

Get the sparse stoichiometry matrix of a model. A feasible solution `x` of a
model `m` is defined as satisfying the equations:

- `stoichiometry(m) * x .== balance(m)`
- `x .>= lbs`
- `y .<= ubs`
- `(lbs, ubs) == bounds(m)
"""
function stoichiometry(a::MetabolicModel)::SparseMat
    _missing_impl_error(stoichiometry, (a,))
end

"""
    bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}

Get the lower and upper solution bounds of a model.
"""
function bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}
    _missing_impl_error(bounds, (a,))
end

"""
    balance(a::MetabolicModel)::SparseVec

Get the sparse balance vector of a model.
"""
function balance(a::MetabolicModel)::SparseVec
    return spzeros(n_metabolites(a))
end

"""
    objective(a::MetabolicModel)::SparseVec

Get the objective vector of the model. Analysis functions, such as
[`flux_balance_analysis`](@ref), are supposed to maximize `dot(objective, x)`
where `x` is a feasible solution of the model.
"""
function objective(a::MetabolicModel)::SparseVec
    _missing_impl_error(objective, (a,))
end


"""
    fluxes(a::MetabolicModel)::Vector{String}

In some models, the [`reactions`](@ref) that correspond to the columns of
[`stoichiometry`](@ref) matrix do not fully represent the semantic contents of
the model; for example, fluxes may be split into forward and reverse reactions,
reactions catalyzed by distinct enzymes, etc. Together with
[`reaction_flux`](@ref) (and [`n_fluxes`](@ref)) this specifies how the
flux is decomposed into individual reactions.

By default (and in most models), fluxes and reactions perfectly correspond.
"""
function fluxes(a::MetabolicModel)::Vector{String}
    reactions(a)
end

function n_fluxes(a::MetabolicModel)::Int
    n_reactions(a)
end

"""
    reaction_flux(a::MetabolicModel)::SparseMat

Retrieve a sparse matrix that describes the correspondence of a solution of the
linear system to the fluxes (see [`fluxes`](@ref) for rationale). Returns a
sparse matrix of size `(n_reactions(a), n_fluxes(a))`. For most models, this is
an identity matrix.
"""
function reaction_flux(a::MetabolicModel)::SparseMat
    nr = n_reactions(a)
    nf = n_fluxes(a)
    nr == nf || _missing_impl_error(reaction_flux, (a,))
    spdiagm(fill(1, nr))
end

"""
    coupling(a::MetabolicModel)::SparseMat

Get a matrix of coupling constraint definitions of a model. By default, there
is no coupling in the models.
"""
function coupling(a::MetabolicModel)::SparseMat
    return spzeros(0, n_reactions(a))
end

"""
    n_coupling_constraints(a::MetabolicModel)::Int

Get the number of coupling constraints in a model.
"""
function n_coupling_constraints(a::MetabolicModel)::Int
    size(coupling(a), 1)
end

"""
    coupling_bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}

Get the lower and upper bounds for each coupling bound in a model, as specified
by `coupling`. By default, the model does not have any coupling bounds.
"""
function coupling_bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}
    return (spzeros(0), spzeros(0))
end

"""
    genes(a::MetabolicModel)::Vector{String}

Return identifiers of all genes contained in the model. By default, there are
no genes.

In SBML, these are usually called "gene products" but we write `genes` for
simplicity.
"""
function genes(a::MetabolicModel)::Vector{String}
    return []
end

"""
    n_genes(a::MetabolicModel)::Int

Return the number of genes in the model (as returned by [`genes`](@ref)). If
you just need the number of the genes, this may be much more efficient than
calling [`genes`](@ref) and measuring the array.
"""
function n_genes(a::MetabolicModel)::Int
    return length(genes(a))
end

"""
    reaction_gene_association(a::MetabolicModel, gene_id::String)::Maybe{GeneAssociation}

Returns the sets of genes that need to be present so that the reaction can work
(technically, a DNF on gene availability, with positive atoms only).

For simplicity, `nothing` may be returned, meaning that the reaction always
takes place. (in DNF, that would be equivalent to returning `[[]]`.)
"""
function reaction_gene_association(
    a::MetabolicModel,
    reaction_id::String,
)::Maybe{GeneAssociation}
    return nothing
end

"""
    reaction_subsystem(model::MetabolicModel, reaction_id::String)::Maybe{String}

Return the subsystem of reaction `reaction_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function reaction_subsystem(model::MetabolicModel, reaction_id::String)::Maybe{String}
    return nothing
end

"""
    reaction_stoichiometry(model::MetaboliteModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid` in the model. The dictionary
maps the metabolite IDs to their stoichiometric coefficients.
"""
function reaction_stoichiometry(m::MetabolicModel, rid::String)::Dict{String,Float64}
    mets = metabolites(m)
    Dict(
        mets[k] => v for
        (k, v) in zip(findnz(stoichiometry(m)[:, first(indexin([rid], reactions(m)))])...)
    )
end

"""
    metabolite_formula(
        a::MetabolicModel,
        metabolite_id::String,
    )::Maybe{MetaboliteFormula}

Return the formula of metabolite `metabolite_id` in `model`.
Return `nothing` in case the formula is not known or irrelevant.
"""
function metabolite_formula(
    model::MetabolicModel,
    metabolite_id::String,
)::Maybe{MetaboliteFormula}
    return nothing
end

"""
    metabolite_charge(model::MetabolicModel, metabolite_id::String)::Maybe{Int}

Return the charge associated with metabolite `metabolite_id` in `model`.
Returns `nothing` if charge not present.
"""
function metabolite_charge(model::MetabolicModel, metabolite_id::String)::Maybe{Int}
    return nothing
end

"""
    metabolite_compartment(model::MetabolicModel, metabolite_id::String)::Maybe{String}

Return the compartment of metabolite `metabolite_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function metabolite_compartment(model::MetabolicModel, metabolite_id::String)::Maybe{String}
    return nothing
end

"""
    reaction_annotations(a::MetabolicModel, reaction_id::String)::Annotations

Return standardized names that may help identifying the reaction. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"Reactome" => ["reactomeID123"]`.
"""
function reaction_annotations(a::MetabolicModel, reaction_id::String)::Annotations
    return Dict()
end

"""
    metabolite_annotations(a::MetabolicModel, metabolite_id::String)::Annotations

Return standardized names that may help to reliably identify the metabolite. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"ChEMBL" => ["123"]` or `"PubChem" => ["CID123", "CID654645645"]`.
"""
function metabolite_annotations(a::MetabolicModel, metabolite_id::String)::Annotations
    return Dict()
end

"""
    gene_annotations(a::MetabolicModel, gene_id::String)::Annotations

Return standardized names that identify the corresponding gene or product. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"PDB" => ["PROT01"]`.
"""
function gene_annotations(a::MetabolicModel, gene_id::String)::Annotations
    return Dict()
end

"""
    reaction_notes(model::MetabolicModel, reaction_id::String)::Notes

Return the notes associated with reaction `reaction_id` in `model`.
"""
function reaction_notes(model::MetabolicModel, reaction_id::String)::Notes
    return Dict()
end

"""
    metabolite_notes(model::MetabolicModel, metabolite_id::String)::Notes

Return the notes associated with metabolite `reaction_id` in `model`.
"""
function metabolite_notes(model::MetabolicModel, metabolite_id::String)::Notes
    return Dict()
end

"""
    gene_notes(model::MetabolicModel, gene_id::String)::Notes

Return the notes associated with the gene `gene_id` in `model`.
"""
function gene_notes(model::MetabolicModel, gene_id::String)::Notes
    return Dict()
end

"""
    enzyme_capacity(model::MetabolicModel)

Return enzyme capacity inequality constraint vector and bound, or nothing
if it doesn't exist in the model.
"""
function enzyme_capacity(model::MetabolicModel)
    #TODO this needs a type
    nothing, nothing
end

"""
    reaction_name(model::MetabolicModel, rid::String)

Return the name of reaction with ID `rid`.
"""
reaction_name(model::MetabolicModel, rid::String) = nothing

"""
    metabolite_name(model::MetabolicModel, mid::String)

Return the name of metabolite with ID `mid`.
"""
metabolite_name(model::MetabolicModel, mid::String) = nothing

"""
    gene_name(model::MetabolicModel, gid::String)

Return the name of gene with ID `gid`.
"""
gene_name(model::MetabolicModel, gid::String) = nothing

"""
    precache!(a::MetabolicModel)::Nothing

Do whatever is feasible to get the model into a state that can be read from
as-quickly-as-possible. This may include e.g. generating helper index
structures and loading delayed parts of the model from disk. The model should
be modified "transparently" in-place. Analysis functions call this right before
applying modifications or converting the model to the optimization model using
[`make_optimization_model`](@ref); usually on the same machine where the
optimizers (and, generally, the core analysis algorithms) will run. The calls
are done in a good hope that the performance will be improved.

By default, it should be safe to do nothing.
"""
function precache!(a::MetabolicModel)::Nothing
    nothing
end
