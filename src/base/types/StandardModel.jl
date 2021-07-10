"""
    mutable struct StandardModel

`StandardModel` is used to store a constraint based metabolic model with meta-information.
Meta-information is defined as annotation details, which include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible, as opposed to `CoreModel` and `CoreModelCoupled`,
which keep the bare neccessities only.
When merging models and keeping meta-information is important, use this as the model type.
If meta-information is not important, use the more efficient core model types.
See [`CoreModel`](@ref) and [`CoreModelCoupled`](@ref) for comparison.

In this model, reactions, metabolites, and genes are stored in ordered dictionaries indexed by each struct's `id` field.
For example, `model.reactions["rxn1_id"]` returns a `Reaction` where the field `id` equals `"rxn1_id"`.
This makes adding and removing reactions efficient.

Note that the stoichiometric matrix (or any other core data, e.g. flux bounds) is not stored directly as in `CoreModel`.
When this model type is used in analysis functions, these core data structures are built from scratch each time an analysis function is called.
This can cause performance issues if you run many small analysis functions sequentially.
Consider using the core model types if performance is critical.

See also: [`Reaction`](@ref), [`Metabolite`](@ref), [`Gene`](@ref)

# Fields
```
id :: String
reactions :: OrderedDict{String, Reaction}
metabolites :: OrderedDict{String, Metabolite}
genes :: OrderedDict{String, Gene}
```

# Example
```
model = load_model(StandardModel, "model_location")
```
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
reactions(model::StandardModel)::StringVecType = collect(keys(model.reactions))

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
metabolites(model::StandardModel)::StringVecType = collect(keys(model.metabolites))

"""
n_metabolites(model::StandardModel)

Return the number of metabolites in `model`.
"""
n_metabolites(model::StandardModel)::Int = length(model.metabolites)

"""
    genes(model::StandardModel)

Return a vector of gene id strings in `model`.
"""
genes(model::StandardModel)::StringVecType = collect(keys(model.genes))

"""
    n_genes(model::StandardModel)

Return the number of genes in `model`.
"""
n_genes(model::StandardModel)::Int = length(model.genes)

"""
    stoichiometry(model::StandardModel)

Return the stoichiometric matrix associated with `model` in sparse format.
"""
function stoichiometry(model::StandardModel)::SparseMat
    S = SparseArrays.spzeros(length(model.metabolites), length(model.reactions))
    met_ids = metabolites(model) # vector of metabolite ids
    rxn_ids = reactions(model)
    for (i, rxn_id) in enumerate(rxn_ids) # column, in order
        for (met_id, coeff) in model.reactions[rxn_id].metabolites
            j = findfirst(==(met_id), met_ids) # row
            if isnothing(j)
                throw(
                    DomainError(
                        met_id,
                        "Metabolite not found in model but occurs in stoichiometry of $(rxn_id). Perhaps it was deleted?",
                    ),
                )
            end
            S[j, i] = coeff
        end
    end
    return S
end

"""
    lower_bounds(model::StandardModel)

Return the lower bounds for all reactions in `model` in sparse format.
"""
function lower_bounds(model::StandardModel)::SparseVec
    sparse([model.reactions[rxn].lb for rxn in reactions(model)])
end

"""
    upper_bounds(model::StandardModel)

Return the upper bounds for all reactions in `model` in sparse format.
Order matches that of the reaction ids returned in `reactions()`.
"""
function upper_bounds(model::StandardModel)::SparseVec
    sparse([model.reactions[rxn].ub for rxn in reactions(model)])
end

"""
    bounds(model::StandardModel)

Return the lower and upper bounds, respectively, for reactions in `model`.
Order matches that of the reaction ids returned in `reactions()`.
"""
function bounds(model::StandardModel)::Tuple{SparseVec,SparseVec}
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
    return sparse([
        model.reactions[rid].objective_coefficient for rid in keys(model.reactions)
    ])
end

"""
    reaction_gene_association(model::StandardModel, id::String)

Return the gene reaction rule in string format for reaction with `id` in `model`.
Return `nothing` if not available.
"""
function reaction_gene_association(model::StandardModel, id::String)::Maybe{GeneAssociation}
    _maybemap(identity, model.reactions[id].grr)
end

"""
    metabolite_formula(model::StandardModel, id::String)

Return the formula of reaction `id` in `model`.
Return `nothing` if not present.
"""
function metabolite_formula(model::StandardModel, id::String)::Maybe{MetaboliteFormula}
    _maybemap(_parse_formula, model.metabolites[id].formula)
end

"""
    metabolite_charge(model::StandardModel, id::String)

Return the charge associated with metabolite `id` in `model`.
Return nothing if not present.
"""
function metabolite_charge(model::StandardModel, id::String)::Maybe{Int}
    model.metabolites[id].charge
end

"""
    metabolite_compartment(model::StandardModel, id::String)

Return compartment associated with metabolite `id` in `model`.
Return `nothing` if not present.
"""
function metabolite_compartment(model::StandardModel, id::String)::Maybe{String}
    model.metabolites[id].compartment
end

"""
    reaction_subsystem(id::String, model::StandardModel)

Return the subsystem associated with reaction `id` in `model`.
Return `nothing` if not present.
"""
function reaction_subsystem(model::StandardModel, id::String)::Maybe{String}
    model.reactions[id].subsystem
end

"""
    metabolite_notes(model::StandardModel, id::String)::Notes

Return the notes associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
function metabolite_notes(model::StandardModel, id::String)::Maybe{Notes}
    model.metabolites[id].notes
end

"""
    metabolite_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
function metabolite_annotations(model::StandardModel, id::String)::Maybe{Annotations}
    model.metabolites[id].annotations
end

"""
    gene_notes(model::StandardModel, id::String)::Notes

Return the notes associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
function gene_notes(model::StandardModel, id::String)::Maybe{Notes}
    model.genes[id].notes
end

"""
    gene_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
function gene_annotations(model::StandardModel, id::String)::Maybe{Annotations}
    model.genes[id].annotations
end

"""
    reaction_notes(model::StandardModel, id::String)::Notes

Return the notes associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
function reaction_notes(model::StandardModel, id::String)::Maybe{Notes}
    model.reactions[id].notes
end

"""
    reaction_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
function reaction_annotations(model::StandardModel, id::String)::Maybe{Annotations}
    model.reactions[id].annotations
end

"""
    reaction_stoichiometry(model::StandardModel, rxn_id::String)::Dict{String, Float64}

Return the reaction equation of reaction with id `rxn_id` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_stoichiometry(m::StandardModel, rxn_id::String)::Dict{String,Float64}
    m.reactions[rxn_id].metabolites
end

"""
Base.convert(::Type{StandardModel}, model::MetabolicModel)

Convert any `MetabolicModel` into a `StandardModel`.
Note, some data loss may occur since only the generic interface is used during
the conversion process.
"""
function Base.convert(::Type{StandardModel}, model::MetabolicModel)
    if typeof(model) == StandardModel
        return model
    end

    id = "" # TODO: add accessor to get model ID
    modelreactions = OrderedDict{String,Reaction}()
    modelmetabolites = OrderedDict{String,Metabolite}()
    modelgenes = OrderedDict{String,Gene}()

    gids = genes(model)
    metids = metabolites(model)
    rxnids = reactions(model)

    for gid in gids
        modelgenes[gid] = Gene(
            gid;
            notes = gene_notes(model, gid),
            annotations = gene_annotations(model, gid),
        ) # TODO: add name accessor
    end

    for mid in metids
        modelmetabolites[mid] = Metabolite(
            mid;
            charge = metabolite_charge(model, mid),
            formula = _maybemap(_unparse_formula, metabolite_formula(model, mid)),
            compartment = metabolite_compartment(model, mid),
            notes = metabolite_notes(model, mid),
            annotations = metabolite_annotations(model, mid),
        )
    end

    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    ocs = objective(model)
    for (i, rid) in enumerate(rxnids)
        rmets = Dict{String,Float64}()
        for (j, stoich) in zip(findnz(S[:, i])...)
            rmets[metids[j]] = stoich
        end
        modelreactions[rid] = Reaction(
            rid;
            metabolites = rmets,
            lb = lbs[i],
            ub = ubs[i],
            grr = reaction_gene_association(model, rid),
            objective_coefficient = ocs[i],
            notes = reaction_notes(model, rid),
            annotations = reaction_annotations(model, rid),
            subsystem = reaction_subsystem(model, rid),
        ) # TODO: add name accessor
    end

    return StandardModel(
        id;
        reactions = modelreactions,
        metabolites = modelmetabolites,
        genes = modelgenes,
    )
end
