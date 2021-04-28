
"""
    abstract type MetabolicModel end

A helper supertype that wraps everything usable as a linear-like model for
COBREXA functions. If you want to use your own type, make it a subtype (so that
the functions typecheck) and add instances for the data accessor methods below.

Required accessor functions are:
1. `reactions`
2. `metabolites`
3. `stoichiometry`
4. `bounds`
5. `objective`
"""
abstract type MetabolicModel end

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
const GeneAssociation = Vector{Vector{String}}
const MetaboliteFormula = Dict{String,Int}
const Annotations = Dict{String,Vector{String}}
const Notes = Dict{String,Vector{String}}
const MatType = AbstractMatrix{Float64}
const VecType = AbstractVector{Float64}
const StringVecType = AbstractVector{String}

_missing_impl_error = (m, a) -> throw(MethodError(m, a))

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

Get the sparse stoichiometry matrix of a model.
"""
function stoichiometry(a::MetabolicModel)::SparseMat
    _missing_impl_error(stoichiometry, (a,))
end

"""
    bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper flux bounds of a model.
"""
function bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}
    _missing_impl_error(bounds, (a,))
end

"""
    balance(a::MetabolicModel)::SparseVec

Get the sparse balance vector of a model (ie. the `b` from `S x = b`).
"""
function balance(a::MetabolicModel)::SparseVec
    return spzeros(n_metabolites(a))
end

"""
    objective(a::MetabolicModel)::SparseVec

Get the objective vector of a model.
"""
function objective(a::MetabolicModel)::SparseVec
    _missing_impl_error(objective, (a,))
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
    coupling_bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper bounds for each coupling bound in a model, as specified
by `coupling`. By default, the model does not have any coupling bounds.
"""
function coupling_bounds(a::MetabolicModel)::Tuple{SparseVec,SparseVec}
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
    metabolite_compartment(model::MetabolicModel, metabolite_id::String)::Maybe{String}

Return the compartment of metabolite `metabolite_id` in `model` if it is assigned. If not, 
return `nothing`. 
"""
function metabolite_compartment(model::MetabolicModel, metabolite_id::String)::Maybe{String}
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
