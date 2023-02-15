"""
$(TYPEDEF)

`ObjectModel` is used to store a constraint based metabolic model with
meta-information.  Meta-information is defined as annotation details, which
include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible, as opposed
to `MatrixModel` and `MatrixModelWithCoupling`, which keep the bare neccessities
only. When merging models and keeping meta-information is important, use this as
the model type.

In this model, reactions, metabolites, and genes are stored in ordered
dictionaries indexed by each struct's `id` field.  For example,
`model.reactions["rxn1_id"]` returns a `Reaction` where the field `id` equals
`"rxn1_id"`.  This makes adding and removing reactions efficient.

Note that the stoichiometric matrix (or any other core data, e.g. flux bounds)
is not stored directly as in [`MatrixModel`](@ref). When this model type is used
in analysis functions, these core data structures are built from scratch each
time an analysis function is called. This can cause performance issues if you
run many small analysis functions sequentially.

See also: [`Reaction`](@ref), [`Metabolite`](@ref), [`Gene`](@ref)

# Example
```
model = load_model(ObjectModel, "my_model.json")
keys(model.reactions)
```

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ObjectModel <: AbstractMetabolicModel
    "Ordered dictionary of reactions."
    reactions::OrderedDict{String,Reaction} = OrderedDict{String,Reaction}()

    "Ordered dictionary of metabolites."
    metabolites::OrderedDict{String,Metabolite} = OrderedDict{String,Metabolite}()

    "Ordered dictionary of genes."
    genes::OrderedDict{String,Gene} = OrderedDict{String,Gene}()

    "Model objective."
    objective::Dict{String,Float64} = Dict{String,Float64}()

    "Machine readable reference to organism embedded via MIRIAM annotation. This should include species name, taxonomy ID, and url to the genome."
    annotations::Annotations = Annotations()

    "Reference information for the model. This should include the DOI and author contact information."
    notes::Notes = Notes()
end

# AbstractMetabolicModel interface follows

Accessors.variables(model::ObjectModel)::StringVecType = collect(keys(model.reactions))

Accessors.n_variables(model::ObjectModel)::Int = length(model.reactions)

Accessors.Internal.@all_variables_are_reactions ObjectModel

Accessors.metabolites(model::ObjectModel)::StringVecType = collect(keys(model.metabolites))

Accessors.n_metabolites(model::ObjectModel)::Int = length(model.metabolites)

Accessors.genes(model::ObjectModel)::StringVecType = collect(keys(model.genes))

Accessors.n_genes(model::ObjectModel)::Int = length(model.genes)

function Accessors.stoichiometry(model::ObjectModel)::SparseMat
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
    rxns = variables(model)
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
    return SparseArrays.sparse(MI, RI, SV, n_metabolites(model), n_variables(model))
end

Accessors.bounds(model::ObjectModel)::Tuple{Vector{Float64},Vector{Float64}} =
    (lower_bounds(model), upper_bounds(model))

Accessors.balance(model::ObjectModel)::SparseVec = spzeros(length(model.metabolites))

Accessors.objective(model::ObjectModel)::SparseVec =
    sparse([get(model.objective, rid, 0.0) for rid in keys(model.reactions)])

function Accessors.reaction_gene_associations(
    model::ObjectModel,
    id::String,
)::Maybe{GeneAssociationsDNF}
    isnothing(model.reactions[id].gene_associations) && return nothing
    [
        collect(keys(rga.gene_product_stoichiometry)) for
        rga in model.reactions[id].gene_associations
    ]
end

Accessors.metabolite_formula(model::ObjectModel, id::String)::Maybe{MetaboliteFormula} =
    maybemap(parse_formula, model.metabolites[id].formula)

Accessors.metabolite_charge(model::ObjectModel, id::String)::Maybe{Int} =
    model.metabolites[id].charge

Accessors.metabolite_compartment(model::ObjectModel, id::String)::Maybe{String} =
    model.metabolites[id].compartment

Accessors.reaction_subsystem(model::ObjectModel, id::String)::Maybe{String} =
    model.reactions[id].subsystem

Accessors.metabolite_notes(model::ObjectModel, id::String)::Maybe{Notes} =
    model.metabolites[id].notes

Accessors.metabolite_annotations(model::ObjectModel, id::String)::Maybe{Annotations} =
    model.metabolites[id].annotations

Accessors.gene_notes(model::ObjectModel, gid::String) = model.genes[gid].notes

Accessors.gene_annotations(model::ObjectModel, id::String)::Maybe{Annotations} =
    model.genes[id].annotations

Accessors.reaction_notes(model::ObjectModel, id::String)::Maybe{Notes} =
    model.reactions[id].notes

Accessors.reaction_annotations(model::ObjectModel, id::String)::Maybe{Annotations} =
    model.reactions[id].annotations

Accessors.reaction_stoichiometry(m::ObjectModel, rid::String)::Dict{String,Float64} =
    m.reactions[rid].metabolites

Accessors.reaction_name(m::ObjectModel, rid::String) = m.reactions[rid].name

Accessors.metabolite_name(m::ObjectModel, mid::String) = m.metabolites[mid].name

Accessors.gene_name(m::ObjectModel, gid::String) = m.genes[gid].name

Accessors.gene_product_molar_mass(model::ObjectModel, gid::String) =
    model.genes[gid].product_molar_mass

Accessors.reaction_isozymes(model::ObjectModel, rid::String) =
    model.reactions[rid].gene_associations

Accessors.gene_product_lower_bound(model::ObjectModel, gid::String) =
    model.genes[gid].product_lower_bound

Accessors.gene_product_upper_bound(model::ObjectModel, gid::String) =
    model.genes[gid].product_upper_bound

Accessors.model_notes(model::ObjectModel)::Notes = model.notes

Accessors.model_annotations(model::ObjectModel)::Annotations = model.annotations

function Base.convert(::Type{ObjectModel}, model::AbstractMetabolicModel)
    if typeof(model) == ObjectModel
        return model
    end

    modelreactions = OrderedDict{String,Reaction}()
    modelmetabolites = OrderedDict{String,Metabolite}()
    modelgenes = OrderedDict{String,Gene}()

    gids = genes(model)
    metids = metabolites(model)
    rxnids = variables(model)

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
            formula = maybemap(unparse_formula, metabolite_formula(model, mid)),
            compartment = metabolite_compartment(model, mid),
            notes = metabolite_notes(model, mid),
            annotations = metabolite_annotations(model, mid),
        )
    end

    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    obj_idxs, obj_vals = findnz(objective(model))
    modelobjective = Dict(k => v for (k, v) in zip(variables(model)[obj_idxs], obj_vals))
    for (i, rid) in enumerate(rxnids)
        rmets = Dict{String,Float64}()
        for (j, stoich) in zip(findnz(S[:, i])...)
            rmets[metids[j]] = stoich
        end
        rgas = reaction_gene_associations(model, rid)
        modelreactions[rid] = Reaction(
            rid;
            name = reaction_name(model, rid),
            metabolites = rmets,
            lower_bound = lbs[i],
            upper_bound = ubs[i],
            gene_associations = isnothing(rgas) ? nothing :
                                [
                Isozyme(; gene_product_stoichiometry = Dict(k => 1.0 for k in rga)) for
                rga in rgas
            ],
            notes = reaction_notes(model, rid),
            annotations = reaction_annotations(model, rid),
            subsystem = reaction_subsystem(model, rid),
        )
    end

    return ObjectModel(;
        reactions = modelreactions,
        metabolites = modelmetabolites,
        genes = modelgenes,
        objective = modelobjective,
        notes = model_notes(model),
        annotations = model_annotations(model),
    )
end
