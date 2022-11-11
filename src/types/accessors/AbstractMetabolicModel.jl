
#
# IMPORTANT
#
# This file provides a list of "officially supported" accessors that should
# work with all subtypes of [`AbstractMetabolicModel`](@ref). Keep this synced with the
# automatically derived methods for [`AbstractModelWrapper`](@ref).
#

#= 
Optimizer interface

These functions interface directly with the optimizer and are primarily used in
make_optimization_model.
=#
"""
$(TYPEDSIGNATURES)

Return `Q`, the matrix used to construct the quadratic component of the
objective in the optimization problem built by
[`make_optimization_model`](@ref). In short, the full objective is `0.5 * x' * Q
* x + q' * x`, where `x` corresponds to the variables in the order returned by
[`variables`](@ref), and `Q` is returned here. Use [`linear_objective`](@ref) to
get `q`.

Returns `nothing` if there is no quadratic component in the model.
"""
function quadratic_objective(a::AbstractMetabolicModel)::Maybe{SparseMat}
    missing_impl_error(quadratic_objective, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `q`, the vector used to construct the linear component of the objective
in the optimization problem built by [`make_optimization_model`](@ref). In
short, the full objective is `0.5 * x' * Q * x + q' * x`, where `x` corresponds
to the variables in the order returned by [`variables`](@ref), and `q` is
returned here. Use [`quadratic_objective`](@ref) to get `Q`.

Returns `nothing` if there is no linear component in the model.
"""
function linear_objective(a::AbstractMetabolicModel)::Maybe{SparseVec}
    missing_impl_error(linear_objective, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `A`, the matrix used to construct the equality constraints in the
optimization problem built by [`make_optimization_model`](@ref). In short the
equality constraints are `A * x == b` where `x` corresponds to the variables in
the order returned by [`variables`](@ref). Use [`balance`](@ref) to get `b`.

Return `nothing` if these constraints are absent in the model.
"""
function stoichiometry(a::AbstractMetabolicModel)::Maybe{SparseMat}
    missing_impl_error(stoichiometry, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `b`, the vector used to construct the equality constraints in the
optimization problem built by [`make_optimization_model`](@ref). In short the
equality constraints are `A * x == b` where `x` corresponds to the variables in
the order returned by [`variables`](@ref). Use [`stoichiometry`](@ref) to get
`A`.

Return `nothing` if these constraints are absent in the model.
"""
function balance(a::AbstractMetabolicModel)::Maybe{SparseVec}
    missing_impl_error(balance, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `C`, the matrix used to construct the inequality constraints in the
optimization problem built by [`make_optimization_model`](@ref). In short, the
inequality constraints are `dₗ ≤ C * x ≤ dᵤ`, where `x` corresponds to the
variables in the order returned by [`variables`](@ref). Use
[`coupling_bounds`](@ref) to get `(dₗ, dᵤ)`.

Return `nothing` if these constraints are absent in the model.
"""
function coupling(a::AbstractMetabolicModel)::Maybe{SparseMat}
    missing_impl_error(coupling, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `(dₗ, dᵤ)`, the bounds used to construct the inequality constraints in
the optimization problem built by [`make_optimization_model`](@ref). In short,
the inequality constraints are `dₗ ≤ C * x ≤ dᵤ`, where `x` corresponds to the
variables in the order returned by [`variables`](@ref). Use [`coupling`](@ref)
to get `C`.

Return `nothing` if these constraints are absent in the model.
"""
function coupling_bounds(a::AbstractMetabolicModel)::Maybe{Tuple{SparseVec, SparseVec}}
    missing_impl_error(coupling_bounds, (a,))
end

"""
$(TYPEDSIGNATURES)

Return `(xₗ, xᵤ)`, the vectors used to construct the simple inequality
constraints (variable bounds) in the optimization problem built by
[`make_optimization_model`](@ref). In short, these simple bounds are `xₗ ≤ x ≤
xᵤ`, where `x` corresponds to the variables in the order returned by
[`variables`](@ref).

Return `nothing` if these constraints are absent in the model.
"""
function bounds(a::AbstractMetabolicModel)::Maybe{Tuple{Vector{Float64},Vector{Float64}}}
    missing_impl_error(bounds, (a,))
end

#= 
The model interface

These functions are used to access the symantic information of the model.
=#

"""
$(TYPEDSIGNATURES)

Return a vector of reaction identifiers in a model. The vector precisely
corresponds to the columns in [`stoichiometry`](@ref) matrix.

For technical reasons, the "reactions" may sometimes not be true reactions but
various virtual and helper pseudo-reactions that are used in the metabolic
modeling, such as metabolite exchanges, separate forward and reverse reactions,
supplies of enzymatic and genetic material and virtual cell volume, etc. To
simplify the view of the model contents use [`reaction_flux`](@ref).
"""
function reactions(a::AbstractMetabolicModel)::Vector{String}
    missing_impl_error(reactions, (a,))
end

"""
$(TYPEDSIGNATURES)

Get the number of reactions in a model.
"""
function n_reactions(a::AbstractMetabolicModel)::Int
    length(reactions(a))
end

"""
$(TYPEDSIGNATURES)

Return a vector of metabolite identifiers in a model. The vector precisely
corresponds to the rows in [`stoichiometry`](@ref) matrix.

As with [`reactions`](@ref)s, some metabolites in models may be virtual,
representing purely technical equality constraints.
"""
function metabolites(a::AbstractMetabolicModel)::Vector{String}
    missing_impl_error(metabolites, (a,))
end

"""
$(TYPEDSIGNATURES)

Get the number of metabolites in a model.
"""
function n_metabolites(a::AbstractMetabolicModel)::Int
    length(metabolites(a))
end

"""
$(TYPEDSIGNATURES)

Return identifiers of all genes contained in the model. By default, there are
no genes.

In SBML, these are usually called "gene products" but we write `genes` for
simplicity.
"""
function genes(a::AbstractMetabolicModel)::Vector{String}
    return String[]
end

"""
$(TYPEDSIGNATURES)

Return the number of genes in the model (as returned by [`genes`](@ref)). If
you just need the number of the genes, this may be much more efficient than
calling [`genes`](@ref) and measuring the array.
"""
function n_genes(a::AbstractMetabolicModel)::Int
    return length(genes(a))
end

"""
$(TYPEDSIGNATURES)

Returns the sets of genes that need to be present so that the reaction can work
(technically, a DNF on gene availability, with positive atoms only).

For simplicity, `nothing` may be returned, meaning that the reaction always
takes place. (in DNF, that would be equivalent to returning `[[]]`.)
"""
function reaction_gene_association(
    a::AbstractMetabolicModel,
    reaction_id::String,
)::Maybe{GeneAssociationsDNF}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the subsystem of reaction `reaction_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function reaction_subsystem(
    model::AbstractMetabolicModel,
    reaction_id::String,
)::Maybe{String}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid` in the model. The dictionary
maps the metabolite IDs to their stoichiometric coefficients.
"""
function reaction_stoichiometry(
    m::AbstractMetabolicModel,
    rid::String,
)::Dict{String,Float64}
    mets = metabolites(m)
    Dict(
        mets[k] => v for
        (k, v) in zip(findnz(stoichiometry(m)[:, first(indexin([rid], reactions(m)))])...)
    )
end

"""
$(TYPEDSIGNATURES)

Return the formula of metabolite `metabolite_id` in `model`.
Return `nothing` in case the formula is not known or irrelevant.
"""
function metabolite_formula(
    model::AbstractMetabolicModel,
    metabolite_id::String,
)::Maybe{MetaboliteFormula}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the charge associated with metabolite `metabolite_id` in `model`.
Returns `nothing` if charge not present.
"""
function metabolite_charge(model::AbstractMetabolicModel, metabolite_id::String)::Maybe{Int}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the compartment of metabolite `metabolite_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function metabolite_compartment(
    model::AbstractMetabolicModel,
    metabolite_id::String,
)::Maybe{String}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return standardized names that may help identifying the reaction. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"Reactome" => ["reactomeID123"]`.
"""
function reaction_annotations(a::AbstractMetabolicModel, reaction_id::String)::Annotations
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return standardized names that may help to reliably identify the metabolite. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"ChEMBL" => ["123"]` or `"PubChem" => ["CID123", "CID654645645"]`.
"""
function metabolite_annotations(
    a::AbstractMetabolicModel,
    metabolite_id::String,
)::Annotations
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return standardized names that identify the corresponding gene or product. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"PDB" => ["PROT01"]`.
"""
function gene_annotations(a::AbstractMetabolicModel, gene_id::String)::Annotations
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return the notes associated with reaction `reaction_id` in `model`.
"""
function reaction_notes(model::AbstractMetabolicModel, reaction_id::String)::Notes
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return the notes associated with metabolite `reaction_id` in `model`.
"""
function metabolite_notes(model::AbstractMetabolicModel, metabolite_id::String)::Notes
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return the notes associated with the gene `gene_id` in `model`.
"""
function gene_notes(model::AbstractMetabolicModel, gene_id::String)::Notes
    return Dict()
end

"""
$(TYPEDSIGNATURES)

Return the name of reaction with ID `rid`.
"""
reaction_name(model::AbstractMetabolicModel, rid::String)::Maybe{String} = nothing

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `mid`.
"""
metabolite_name(model::AbstractMetabolicModel, mid::String)::Maybe{String} = nothing

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `gid`.
"""
gene_name(model::AbstractMetabolicModel, gid::String)::Maybe{String} = nothing

#= 
Output interface

These functions are used to correctly associate optimizer variables with
semantic model information.
=#
"""
$(TYPEDSIGNATURES)

In some models, the [`reactions`](@ref) that correspond to the columns of
[`stoichiometry`](@ref) matrix do not fully represent the semantic contents of
the model; for example, fluxes may be split into forward and reverse reactions,
reactions catalyzed by distinct enzymes, etc. By using the semantic variable
types, [`EnzymeAbundances`](@ref), [`ReactionFluxes`](@ref), etc. with
[`map_optimizer_to_semantic_variables`](@ref) values of the optimization model
can be mapped to meaningful quantities.

# Usage
```
gm = GeckoModel(...)
opt_model = flux_balance_analysis(gm, ...)

enzyme_abundances = Dict(
    semantic_variables(gm, EnzymeAbundances) .=> map_optimizer_to_semantic_variables(gm, EnzymeAbundances, value.(opt_model))
)
```
"""
function semantic_variables(a::AbstractMetabolicModel, b::AbstractSemanticVariables)::Vector{String}
    reactions(a)
end

function n_semantic_variables(a::AbstractMetabolicModel, b::AbstractSemanticVariables)::Int
    n_reactions(a)
end

"""
$(TYPEDSIGNATURES)

The ids of the variables solved for by the optimization model. This corresponds
to the ids of the columns returned by [`stoichiometry`](@ref). These do not
necessarily need to correspond to semantically meaningful variables, but often
they do.
"""
function optimizer_variables(a::AbastractMetabolicModel)
    reactions(a)
end

function n_optimizer_variables(a::AbstractMetabolicModel)
    n_reactions(a)
end

"""
$(TYPEDSIGNATURES)

Transforms the output of the optimizer, `opt_output`, to some semantically
meaningful variables, `b`, for model, `a`. Typically this will just be the
identity transform. 

This function becomes useful when special transformations on variables need to
happen, e.g. for enzyme constrained models, bidirectional reactions are
typically converted into two unidirectional reactions. Consequently, there are
two optimizer variables, but still only one semantically important variable.
"""
function map_optimizer_to_semantic_variables(
    a::AbstractMetabolicModel, 
    b::AbstractSemanticVariables, 
    opt_output::Vector{Float64},
)::Vector{Float64}
    opt_output
end

#= 
Miscellaneous interface.

These functions perform various other useful tasks.
=#
"""
$(TYPEDSIGNATURES)

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
function precache!(a::AbstractMetabolicModel)::Nothing
    nothing
end
