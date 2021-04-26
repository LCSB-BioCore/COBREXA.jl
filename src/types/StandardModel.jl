"""
    mutable struct StandardModel

`StandardModel` is used to store a constraint based metabolic model with meta-information.
Meta-information is defined as annotation details, which include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible, cf. `CoreModel` and `CoreModelCoupled`.
When merging models and keeping meta-information is important, use this as the model type. 
If meta-information is not important, use the more efficient core model types. 
See [`CoreModel`](@ref) and [`CoreModelCoupled`](@ref) for comparison.

In this model, reactions, metabolites, and genes are stored in dictionaries indexed by each struct's `id` field.
For example, `model.reactions["rxn1_id"]` returns a `Reaction` that is indexed by the field `id` which equals `"rxn1_id"`.
This makes adding and removing reactions efficient.   

However, note that the stoichiometric matrix (or any other core data, e.g. flux bounds) is not stored directly as in `CoreModel`. 
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
        id = "";
        reactions = OrderedDict{String,Reaction}(),
        metabolites = OrderedDict{String,Metabolite}(),
        genes = OrderedDict{String,Gene}(),
    ) = new(id, reactions, metabolites, genes)
end

# MetabolicModel interface follows
"""
reactions(model::StandardModel)

Return a vector of reaction id strings contained in `model`.
The order of reaction ids returned here matches the order used to construct the
stoichiometric matrix.
"""
reactions(model::StandardModel)::Vector{String} = [r_id for r_id in keys(model.reactions)]

"""
    n_reactions(model::StandardModel)

Return the number of reactions contained in `model`.
"""
n_reactions(model::StandardModel)::Int = length(model.reactions)

"""
    metabolites(model::StandardModel)

Return a vector of metabolite id strings contained in `model`.
The order of metabolite strings returned here matches the order used to construct
the stoichiometric matrix.
"""
metabolites(model::StandardModel)::Vector{String} =
    [m_id for m_id in keys(model.metabolites)]

"""
n_metabolites(model::StandardModel)

Return the number of metabolites in `model`.
"""
n_metabolites(model::StandardModel)::Int = length(model.metabolites)

"""
    genes(model::StandardModel)

Return a vector of gene id strings in `model`.
"""
genes(model::StandardModel)::Vector{String} = [g_id for g_id in keys(model.genes)]

"""
    n_genes(model::StandardModel)

Return the number of genes in `model`.
"""
n_genes(model::StandardModel)::Int = length(model.genes)

"""
    stoichiometry(model::StandardModel)

Return the stoichiometric matrix associated with `model` in sparse format.
"""
function stoichiometry(model::StandardModel)
    S = SparseArrays.spzeros(length(model.metabolites), length(model.reactions))
    met_ids = metabolites(model) # vector of metabolite ids
    rxn_ids = reactions(model)
    for (i, rxn_id) in enumerate(rxn_ids) # column, in order
        for (met_id, coeff) in model.reactions[rxn_id].metabolites
            j = findfirst(x -> x == met_id, met_ids) # row
            S[j, i] = coeff
        end
    end
    return S
end

"""
    lower_bounds(model::StandardModel)

Return the lower bounds for all reactions in `model` in sparse format.
"""
function lower_bounds(model::StandardModel)
    sparse([model.reactions[rxn].lb for rxn in reactions(model)])
end

"""
    upper_bounds(model::StandardModel)

Return the upper bounds for all reactions in `model` in sparse format.
"""
function upper_bounds(model::StandardModel)
    sparse([model.reactions[rxn].ub for rxn in reactions(model)])
end

"""
    bounds(model::StandardModel)

Return the lower and upper bounds, respectively, for reactions in `model`.
"""
function bounds(model::StandardModel)
    ubs = upper_bounds(model)
    lbs = lower_bounds(model)
    return lbs, ubs
end

"""
    balance(model::StandardModel)

Return the balance of the linear problem, i.e. b in Sv = 0 where S is the stoichiometric matrix
and v is the flux vector.
"""
balance(model::StandardModel)::SparseVec = spzeros(length(model.metabolites))

"""
    objective(model::StandardModel)

Return sparse objective vector for `model`.
"""
function objective(model::StandardModel)::SparseVec
    obj_arr = SparseArrays.spzeros(length(model.reactions))
    for (i, r_id) in enumerate(reactions(model))
        if model.reactions[r_id].objective_coefficient != 0.0
            obj_arr[i] = 1.0 # could have multiple objective coefficients
        end
    end
    nnz(obj_arr) == 0 && (@warn "No objective found.")
    return obj_arr
end

"""
    reaction_gene_association(id::String, model::StandardModel)

Return the gene reaction rule for reaction with `id` in `model`.
"""
function reaction_gene_association(id::String, model::StandardModel)
    maybemap(_unparse_grr, model.reactions[id].grr) # this needs to be mapped to a boolean but that isn't implemented yet
end

"""
    formula(id::String, model::StandardModel)

Return the formula of reaction `id` in `model`.
"""
function formula(id::String, model::StandardModel)
    maybemap(get_atoms, model.metabolites[id].formula)
end

"""
    charge(id::String, model::StandardModel)

Return the charge associated with metabolite `id` in `model`.
"""
function charges(id::String, model::StandardModel)
    model.metabolites[id].charge
end

"""
    metabolite_chemistry(model::StandardModel)

Return the formula and charge associated with metabolite `id` in `model`.
"""
function metabolite_chemistry(id::String, model::StandardModel)
    (formula(id, model), charge(id, model))
end

"""
    metabolite_compartment(id::String, model::StandardModel)

Return compartment associated with metabolite `id` in `model`.
"""
function metabolite_compartment(id::String, model::StandardModel)
    model.metabolites[id].compartment
end

"""
    reaction_subsystem(id::String, model::StandardModel)

Return the subsystem associated with reaction `id` in `model`.
"""
function reaction_subsystem(id::String, model::StandardModel)
    model.reactions[id].subsystem
end

"""
    metabolite_notes(id::String, model::StandardModel)

Return the notes associated with metabolite `id` in `model`.
"""
function metabolite_notes(id::String, model::StandardModel)
    model.metabolites[id].notes
end

"""
    metabolite_annotation(id::String, model::StandardModel)

Return the annotation associated with metabolite `id` in `model`.
"""
function metabolite_annotation(id::String, model::StandardModel)
    model.metabolites[id].annotation
end

"""
    gene_notes(id::String, model::StandardModel)

Return the notes associated with gene `id` in `model`.
"""
function gene_notes(id::String, model::StandardModel)
    [g.notes for g in model.genes]
end

"""
    gene_annotations(id::String, model::StandardModel)

Return the annotation associated with gene `id` in `model`.
"""
function gene_annotations(id::String, model::StandardModel)
    [g.annotation for g in model.genes]
end

"""
    reaction_notes(id::String, model::StandardModel)

Return the notes associated with reaction `id` in `model`.
"""
function reaction_notes(id::String, model::StandardModel)
    [r.notes for r in model.reactions]
end

"""
    reaction_annotations(id::String, model::StandardModel)

Return the annotation associated with reaction `id` in `model`.
"""
function reaction_annotations(id::String, model::StandardModel)
    [r.annotation for r in model.reactions]
end

function Base.convert(::Type{StandardModel}, model::MetabolicModel)
    id = "" # model_id(model), add accessor
    modelreactions = OrderedDict{String,Reaction}()
    modelmetabolites = OrderedDict{String,Metabolite}()
    modelgenes = OrderedDict{String,Gene}()

    gids = genes(model)
    metids = metabolites(model)
    rxnids = reactions(model)

    for gid in gids
        modelgenes[gid] = Gene(gid;) # NB: add more accessors
    end

    for mid in metids
        f, c = metabolite_chemistry(model, mid)
        fstr = join([k * string(v) for (k, v) in f])
        modelmetabolites[mid] = Metabolite(mid; charge=c, formula=fstr)
    end

    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    for (i, rid) in enumerate(rxnids)
        rmets = Dict{String,Float64}()
        for (j, stoich) in zip(findnz(S[:, i])...)
            rmets[metids[j]] = stoich
        end
        modelreactions[rid] = Reaction(rid; metabolites = rmets, lb=lbs[i], ub=ubs[i], grr=reaction_gene_association(model, rid)) # NB: add more accessors
    end

    return StandardModel(
        id;
        reactions = modelreactions,
        metabolites = modelmetabolites,
        genes = modelgenes,
    )
end
