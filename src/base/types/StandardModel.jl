"""
    mutable struct StandardModel

`StandardModel` is used to store a constraint based metabolic model with
meta-information.  Meta-information is defined as annotation details, which
include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible, as opposed
to `CoreModel` and `CoreModelCoupled`, which keep the bare neccessities only.
When merging models and keeping meta-information is important, use this as the
model type.  If meta-information is not important, use the more efficient core
model types.  See [`CoreModel`](@ref) and [`CoreModelCoupled`](@ref) for
comparison.

In this model, reactions, metabolites, and genes are stored in ordered
dictionaries indexed by each struct's `id` field.  For example,
`model.reactions["rxn1_id"]` returns a `Reaction` where the field `id` equals
`"rxn1_id"`.  This makes adding and removing reactions efficient.

Note that the stoichiometric matrix (or any other core data, e.g. flux bounds)
is not stored directly as in `CoreModel`.  When this model type is used in
analysis functions, these core data structures are built from scratch each time
an analysis function is called.  This can cause performance issues if you run
many small analysis functions sequentially.  Consider using the core model
types if performance is critical.

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
model = load_model(StandardModel, "my_model.json")
keys(model.reactions)
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
    n_entries = 0
    for (_, r) in model.reactions
        for _ in r.metabolites
            n_entries += 1
        end
    end

    MI = Vector{Int}()
    RI = Vector{Int}()
    SV = Vector{Float64}()
    sizehint!(MI, n_entries)
    sizehint!(RI, n_entries)
    sizehint!(SV, n_entries)

    # establish the ordering
    rxns = reactions(model)
    met_idx = Dict(mid => i for (i, mid) in enumerate(metabolites(model)))

    # fill the matrix entries
    for (ridx, rid) in enumerate(rxns)
        for (mid, coeff) in model.reactions[rid].metabolites
            haskey(met_idx, mid) || throw(
                DomainError(
                    mid,
                    "Metabolite $(mid) not found in model but occurs in stoichiometry of $(rid)",
                ),
            )
            push!(MI, met_idx[mid])
            push!(RI, ridx)
            push!(SV, coeff)
        end
    end
    return SparseArrays.sparse(MI, RI, SV, n_metabolites(model), n_reactions(model))
end

"""
    lower_bounds(model::StandardModel)::Vector{Float64}

Return the lower bounds for all reactions in `model` in sparse format.
"""
lower_bounds(model::StandardModel)::Vector{Float64} =
    sparse([model.reactions[rxn].lb for rxn in reactions(model)])

"""
    upper_bounds(model::StandardModel)::Vector{Float64}

Return the upper bounds for all reactions in `model` in sparse format.
Order matches that of the reaction ids returned in `reactions()`.
"""
upper_bounds(model::StandardModel)::Vector{Float64} =
    sparse([model.reactions[rxn].ub for rxn in reactions(model)])

"""
    bounds(model::StandardModel)::Tuple{Vector{Float64},Vector{Float64}}

Return the lower and upper bounds, respectively, for reactions in `model`.
Order matches that of the reaction ids returned in `reactions()`.
"""
bounds(model::StandardModel)::Tuple{Vector{Float64},Vector{Float64}} =
    (lower_bounds(model), upper_bounds(model))

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
objective(model::StandardModel)::SparseVec =
    sparse([model.reactions[rid].objective_coefficient for rid in keys(model.reactions)])

"""
    reaction_gene_association(model::StandardModel, id::String)

Return the gene reaction rule in string format for reaction with `id` in `model`.
Return `nothing` if not available.
"""
reaction_gene_association(model::StandardModel, id::String)::Maybe{GeneAssociation} =
    _maybemap(identity, model.reactions[id].grr)

"""
    metabolite_formula(model::StandardModel, id::String)

Return the formula of reaction `id` in `model`.
Return `nothing` if not present.
"""
metabolite_formula(model::StandardModel, id::String)::Maybe{MetaboliteFormula} =
    _maybemap(_parse_formula, model.metabolites[id].formula)

"""
    metabolite_charge(model::StandardModel, id::String)

Return the charge associated with metabolite `id` in `model`.
Return nothing if not present.
"""
metabolite_charge(model::StandardModel, id::String)::Maybe{Int} =
    model.metabolites[id].charge

"""
    metabolite_compartment(model::StandardModel, id::String)

Return compartment associated with metabolite `id` in `model`.
Return `nothing` if not present.
"""
metabolite_compartment(model::StandardModel, id::String)::Maybe{String} =
    model.metabolites[id].compartment

"""
    reaction_subsystem(id::String, model::StandardModel)

Return the subsystem associated with reaction `id` in `model`.
Return `nothing` if not present.
"""
reaction_subsystem(model::StandardModel, id::String)::Maybe{String} =
    model.reactions[id].subsystem

"""
    metabolite_notes(model::StandardModel, id::String)::Notes

Return the notes associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
metabolite_notes(model::StandardModel, id::String)::Maybe{Notes} =
    model.metabolites[id].notes

"""
    metabolite_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
metabolite_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.metabolites[id].annotations

"""
    gene_notes(model::StandardModel, id::String)::Notes

Return the notes associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
gene_notes(model::StandardModel, id::String)::Maybe{Notes} = model.genes[id].notes

"""
    gene_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
gene_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.genes[id].annotations

"""
    reaction_notes(model::StandardModel, id::String)::Notes

Return the notes associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
reaction_notes(model::StandardModel, id::String)::Maybe{Notes} = model.reactions[id].notes

"""
    reaction_annotations(model::StandardModel, id::String)::Annotations

Return the annotation associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
reaction_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.reactions[id].annotations

"""
    reaction_stoichiometry(model::StandardModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::StandardModel, rid::String)::Dict{String,Float64} =
    m.reactions[rid].metabolites

"""
    reaction_name(m::StandardModel, rid::String)

Return the name of reaction with ID `id`.
"""
reaction_name(m::StandardModel, rid::String) = m.reactions[rid].name

"""
    metabolite_name(m::StandardModel, mid::String)

Return the name of metabolite with ID `id`.
"""
metabolite_name(m::StandardModel, mid::String) = m.metabolites[mid].name

"""
    gene_name(m::StandardModel, gid::String)

Return the name of gene with ID `id`.
"""
gene_name(m::StandardModel, gid::String) = m.genes[gid].name

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
            name = gene_name(model, gid),
            notes = gene_notes(model, gid),
            annotations = gene_annotations(model, gid),
        )
    end

    for mid in metids
        modelmetabolites[mid] = Metabolite(
            mid;
            name = metabolite_name(model, mid),
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
            name = reaction_name(model, rid),
            metabolites = rmets,
            lb = lbs[i],
            ub = ubs[i],
            grr = reaction_gene_association(model, rid),
            objective_coefficient = ocs[i],
            notes = reaction_notes(model, rid),
            annotations = reaction_annotations(model, rid),
            subsystem = reaction_subsystem(model, rid),
        )
    end

    return StandardModel(
        id;
        reactions = modelreactions,
        metabolites = modelmetabolites,
        genes = modelgenes,
    )
end

#TODO generalize these to other model types

"""
    reaction_bounds(model::StandardModel, rid::String)

Return lower and upper bounds for `rid` in `model`.
"""
function reaction_bounds(model::StandardModel, rid::String)
    model.reactions[rid].lb, model.reactions[rid].ub
end

"""
    is_reaction_reversible(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is reversible.
"""
function is_reaction_reversible(model::StandardModel, rid::String)
    lb, ub = reaction_bounds(model, rid)
    lb < 0 && ub > 0
end

"""
    is_reaction_forward_only(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is forward only.
"""
function is_reaction_forward_only(model::StandardModel, rid::String)
    lb, ub = reaction_bounds(model, rid)
    lb >= 0 && ub > 0
end

"""
    is_reaction_backward_only(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is backward only.
"""
function is_reaction_backward_only(model::StandardModel, rid::String)
    lb, ub = reaction_bounds(model, rid)
    lb < 0 && ub <= 0
end

"""
    is_reaction_unidirectional(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is unidirectional.
"""
function is_reaction_unidirectional(model::StandardModel, rid::String)
    is_reaction_forward_only(model, rid) || is_reaction_backward_only(model, rid)
end

"""
    is_reaction_blocked(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is blocked.
"""
function is_reaction_blocked(model::StandardModel, rid::String)
    lb, ub = reaction_bounds(model, rid)
    lb == ub == 0
end

"""
    has_reaction_isozymes(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is catalyzed by multiple enzymes,
i.e. it has isozymes according to the gene reaction rules.
"""
function has_reaction_isozymes(model::StandardModel, rid::String)
    length(reaction_gene_association(model, rid)) > 1
end

"""
    reaction_has_grr(model::StandardModel, rid::String)

Check if reaction `rid` in `model` has a gene reaction rule entry.
"""
function has_reaction_grr(model::StandardModel, rid::String)
    #TODO simplify this once COBREXA enforces universal rules for GRR representation
    !isnothing(reaction_gene_association(model, rid)) &&
        reaction_gene_association(model, rid) != [[]] &&
        !isempty(first(reaction_gene_association(model, rid)))
end
