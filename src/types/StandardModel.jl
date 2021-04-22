"""
`StandardModel` is used to store a constraint based metabolic model with meta-information.
Meta-information is defined as annotation details, which include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible.
When merging models and keeping meta-information is important, use this as the output type. 
If meta-information is not important, use the more efficient core model types. 
See [`CoreModel`](@ref) and [`CoreModelCoupled`](@ref) for comparison.

In this model, reactions, metabolites, and genes are stored in dictionaries indexed by each structs `id` field.
For example, `model.reactions["rxn1_id"]` returns a `Reaction` with index `rxn1_id`.
This makes adding and removing reactions efficient.   

However, note that the stoichiometric matrix (or any other core data, e.g. flux bounds) is not stored directly. 
When this model type is used in analysis functions, these core data structures are built from scratch.
This can cause performance issues if you run many small analysis functions sequentially. 
Consider using the core model types if performance is critical.

See also: [`Reaction`](@ref), [`Metabolite`](@ref), [`Gene`](@ref)

# Fields
````
id :: String
reactions :: OrderedDict{String, Reaction}
metabolites :: OrderedDict{String, Metabolite}
genes :: OrderedDict{String, Gene}
````
"""
mutable struct StandardModel <: MetabolicModel
    id::String
    reactions::OrderedDict{String,Reaction}
    metabolites::OrderedDict{String,Metabolite}
    genes::OrderedDict{String,Gene}

    StandardModel(
        id = "",
        reactions = OrderedDict{String,Reaction}(),
        metabolites = OrderedDict{String,Metabolite}(),
        genes = OrderedDict{String,Gene}(),
    ) = new(id, reactions, metabolites, genes)
end

# MetabolicModel interface follows
reactions(model::StandardModel)::Vector{String} = [r_id for r_id in keys(model.reactions)]
n_reactions(model::StandardModel)::Int = length(model.reactions)

metabolites(model::StandardModel)::Vector{String} = [m_id for m_id in keys(model.metabolites)]
n_metabolites(model::StandardModel)::Int = length(model.metabolites)

genes(model::StandardModel)::Vector{String} = [g_id for g_id in keys(model.genes)]
n_genes(model::StandardModel)::Int = length(model.genes)

function stoichiometry(model::StandardModel)
    S = SparseArrays.spzeros(length(model.metabolites), length(model.reactions))
    met_ids = metabolites(model) # vector of metabolite ids
    for (i, rxn_id) in enumerate(reactions(model)) # column, in order
        for (met_id, coeff) in model.reactions[rxn_id].metabolites
            j = findfirst(x -> x == met_id, met_ids) # row
            isnothing(j) ?
            (@error "S matrix construction error: $(met_id) not defined."; return nothing) :
            nothing
            S[j, i] = coeff
        end
    end
    return S
end

function lower_bounds(model::StandardModel)
    [model.reactions[rxn].lb for rxn in reactions(model)]
end

function upper_bounds(model::StandardModel)
    [model.reactions[rxn].ub for rxn in reactions(model)]
end

function bounds(model::StandardModel)
    ubs = upper_bounds(model)
    lbs = lower_bounds(model)
    return lbs, ubs
end

balance(model::StandardModel)::SparseVec = spzeros(length(model.metabolites))

function objective(model::StandardModel)::SparseVec
    obj_arr = SparseArrays.spzeros(length(model.reactions))
    j = -1
    for (i, r) in enumerate(reactions(model))
        if model[rxn].objective_coefficient != 0.0
            j = i
            obj_arr[j] = 1.0 # could have multiple objective coefficients
            break
        end
    end
    if j != -1 # objective assigned, otherwise return array of 0s
        obj_arr[j] = 1.0
    else
        nnz(c) == 0 && (@warn "No objective found.")
    end
    return obj_arr
end

function gene_associations(model::StandardModel)
    grrs = [_unparse_grr(r.grr) for r in values(model.reactions)]
end

function formulas(model::StandardModel)
    [_formula_to_atoms(m.formula) for m in model.metabolites]
end

function charges(model::StandardModel)
    [m.charge for m in model.metabolites]
end

function metabolite_chemistry(model::StandardModel)
    (formulas(model), charges(model))
end

function metabolite_compartments(model::StandardModel)
    [m.compartment for m in model.metabolites]
end

function reaction_subsystems(model::StandardModel)
    [r.subsystem for r in model.reactions]
end

function metabolite_notes(model::JSONModel)
    [m.notes for m in model.metabolites]
end

function metabolite_annotations(model::JSONModel)
    [m.annotation for m in model.metabolites]
end

function gene_notes(model::JSONModel)
    [g.notes for g in model.genes]
end

function gene_annotations(model::JSONModel)
    [g.annotation for g in model.genes]
end

function reaction_notes(model::JSONModel)
    [r.notes for r in model.reactions]
end

function reaction_annotations(model::JSONModel)
    [r.annotation for r in model.reactions]
end
