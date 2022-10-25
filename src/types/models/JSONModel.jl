"""
$(TYPEDEF)

A struct used to store the contents of a JSON model, i.e. a model read from a
file ending with `.json`. These model files typically store all the model data
in arrays of JSON objects (represented in Julia as vectors of dictionaries).

Usually, not all of the fields of the input JSON can be easily represented when
converting to other models, care should be taken to avoid losing information.

Direct work with the `json` structure is not very efficient; the model
structure therefore caches some of the internal structure in the extra fields.
The single-parameter [`JSONModel`](@ref) constructor creates these caches
correctly from the `json`. The model structure is designed as read-only, and
changes in `json` invalidate the cache.

# Example
````
model = load_json_model("some_model.json")
model.json # see the actual underlying JSON
reactions(model) # see the list of reactions
````

# Fields
$(TYPEDFIELDS)
"""
struct JSONModel <: AbstractMetabolicModel
    json::Dict{String,Any}
    rxn_index::Dict{String,Int}
    rxns::Vector{Any}
    met_index::Dict{String,Int}
    mets::Vector{Any}
    gene_index::Dict{String,Int}
    genes::Vector{Any}
end

_json_rxn_name(r, i) = string(get(r, "id", "rxn$i"))
_json_met_name(m, i) = string(get(m, "id", "met$i"))
_json_gene_name(g, i) = string(get(g, "id", "gene$i"))

JSONModel(json::Dict{String,Any}) = begin
    rkey = guesskey(keys(json), constants.keynames.rxns)
    isnothing(rkey) && throw(DomainError(keys(json), "JSON model has no reaction keys"))
    rs = json[rkey]

    mkey = guesskey(keys(json), constants.keynames.mets)
    ms = json[mkey]
    isnothing(mkey) && throw(DomainError(keys(json), "JSON model has no metabolite keys"))

    gkey = guesskey(keys(json), constants.keynames.genes)
    gs = isnothing(gkey) ? [] : json[gkey]

    JSONModel(
        json,
        Dict(_json_rxn_name(r, i) => i for (i, r) in enumerate(rs)),
        rs,
        Dict(_json_met_name(m, i) => i for (i, m) in enumerate(ms)),
        ms,
        Dict(_json_gene_name(g, i) => i for (i, g) in enumerate(gs)),
        gs,
    )
end

function _parse_annotations(x)::Annotations
    Annotations([k => if typeof(v) == String
        [v]
    else
        convert(Vector{String}, v)
    end for (k, v) in x])
end

_parse_notes(x)::Notes = _parse_annotations(x)

Accessors.n_reactions(model::JSONModel) = length(model.rxns)
Accessors.n_metabolites(model::JSONModel) = length(model.mets)
Accessors.n_genes(model::JSONModel) = length(model.genes)

"""
$(TYPEDSIGNATURES)

Extract reaction names (stored as `.id`) from JSON model.
"""
Accessors.reactions(model::JSONModel) =
    [_json_rxn_name(r, i) for (i, r) in enumerate(model.rxns)]

"""
$(TYPEDSIGNATURES)

Extract metabolite names (stored as `.id`) from JSON model.
"""
Accessors.metabolites(model::JSONModel) =
    [_json_met_name(m, i) for (i, m) in enumerate(model.mets)]

"""
$(TYPEDSIGNATURES)

Extract gene names from a JSON model.
"""
Accessors.genes(model::JSONModel) =
    [_json_gene_name(g, i) for (i, g) in enumerate(model.genes)]

"""
$(TYPEDSIGNATURES)

Get the stoichiometry. Assuming the information is stored in reaction object
under key `.metabolites`.
"""
function Accessors.stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    met_ids = metabolites(model)

    n_entries = 0
    for r in model.rxns
        for _ in r["metabolites"]
            n_entries += 1
        end
    end

    MI = Vector{Int}()
    RI = Vector{Int}()
    SV = Vector{Float64}()
    sizehint!(MI, n_entries)
    sizehint!(RI, n_entries)
    sizehint!(SV, n_entries)

    for (i, rid) in enumerate(rxn_ids)
        r = model.rxns[model.rxn_index[rid]]
        for (mid, coeff) in r["metabolites"]
            haskey(model.met_index, mid) || throw(
                DomainError(
                    met_id,
                    "Unknown metabolite found in stoichiometry of $(rxn_ids[i])",
                ),
            )

            push!(MI, model.met_index[mid])
            push!(RI, i)
            push!(SV, coeff)
        end
    end
    return SparseArrays.sparse(MI, RI, SV, length(met_ids), length(rxn_ids))
end

"""
$(TYPEDSIGNATURES)

Get the bounds for reactions, assuming the information is stored in
`.lower_bound` and `.upper_bound`.
"""
Accessors.bounds(model::JSONModel) = (
    [get(rxn, "lower_bound", -constants.default_reaction_bound) for rxn in model.rxns],
    [get(rxn, "upper_bound", constants.default_reaction_bound) for rxn in model.rxns],
)

"""
$(TYPEDSIGNATURES)

Collect `.objective_coefficient` keys from model reactions.
"""
Accessors.objective(model::JSONModel) =
    sparse([float(get(rxn, "objective_coefficient", 0.0)) for rxn in model.rxns])

"""
$(TYPEDSIGNATURES)

Parses the `.gene_reaction_rule` from reactions.
"""
Accessors.reaction_gene_association(model::JSONModel, rid::String) = maybemap(
    parse_grr,
    get(model.rxns[model.rxn_index[rid]], "gene_reaction_rule", nothing),
)

"""
$(TYPEDSIGNATURES)

Parses the `.subsystem` out from reactions.
"""
Accessors.reaction_subsystem(model::JSONModel, rid::String) =
    get(model.rxns[model.rxn_index[rid]], "subsystem", nothing)

"""
$(TYPEDSIGNATURES)

Parse and return the metabolite `.formula`
"""
Accessors.metabolite_formula(model::JSONModel, mid::String) =
    maybemap(parse_formula, get(model.mets[model.met_index[mid]], "formula", nothing))

"""
$(TYPEDSIGNATURES)

Return the metabolite `.charge`
"""
Accessors.metabolite_charge(model::JSONModel, mid::String) =
    get(model.mets[model.met_index[mid]], "charge", 0)

"""
$(TYPEDSIGNATURES)

Return the metabolite `.compartment`
"""
Accessors.metabolite_compartment(model::JSONModel, mid::String) =
    get(model.mets[model.met_index[mid]], "compartment", nothing)

"""
$(TYPEDSIGNATURES)

Gene annotations from the [`JSONModel`](@ref).
"""
Accessors.gene_annotations(model::JSONModel, gid::String)::Annotations = maybemap(
    _parse_annotations,
    get(model.genes[model.gene_index[gid]], "annotation", nothing),
)

"""
$(TYPEDSIGNATURES)

Gene notes from the [`JSONModel`](@ref).
"""
Accessors.gene_notes(model::JSONModel, gid::String)::Notes =
    maybemap(_parse_notes, get(model.genes[model.gene_index[gid]], "notes", nothing))

"""
$(TYPEDSIGNATURES)

Reaction annotations from the [`JSONModel`](@ref).
"""
Accessors.reaction_annotations(model::JSONModel, rid::String)::Annotations = maybemap(
    _parse_annotations,
    get(model.rxns[model.rxn_index[rid]], "annotation", nothing),
)

"""
$(TYPEDSIGNATURES)

Reaction notes from the [`JSONModel`](@ref).
"""
Accessors.reaction_notes(model::JSONModel, rid::String)::Notes =
    maybemap(_parse_notes, get(model.rxns[model.rxn_index[rid]], "notes", nothing))

"""
$(TYPEDSIGNATURES)

Metabolite annotations from the [`JSONModel`](@ref).
"""
Accessors.metabolite_annotations(model::JSONModel, mid::String)::Annotations = maybemap(
    _parse_annotations,
    get(model.mets[model.met_index[mid]], "annotation", nothing),
)

"""
$(TYPEDSIGNATURES)

Metabolite notes from the [`JSONModel`](@ref).
"""
Accessors.metabolite_notes(model::JSONModel, mid::String)::Notes =
    maybemap(_parse_notes, get(model.mets[model.met_index[mid]], "notes", nothing))

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
Accessors.reaction_stoichiometry(model::JSONModel, rid::String)::Dict{String,Float64} =
    model.rxns[model.rxn_index[rid]]["metabolites"]

"""
$(TYPEDSIGNATURES)

Return the name of reaction with ID `rid`.
"""
Accessors.reaction_name(model::JSONModel, rid::String) =
    get(model.rxns[model.rxn_index[rid]], "name", nothing)

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `mid`.
"""
Accessors.metabolite_name(model::JSONModel, mid::String) =
    get(model.mets[model.met_index[mid]], "name", nothing)

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `gid`.
"""
Accessors.gene_name(model::JSONModel, gid::String) =
    get(model.genes[model.gene_index[gid]], "name", nothing)

"""
$(TYPEDSIGNATURES)

Convert any [`AbstractMetabolicModel`](@ref) to [`JSONModel`](@ref).
"""
function Base.convert(::Type{JSONModel}, mm::AbstractMetabolicModel)
    if typeof(mm) == JSONModel
        return mm
    end

    rxn_ids = reactions(mm)
    met_ids = metabolites(mm)
    gene_ids = genes(mm)
    S = stoichiometry(mm)
    lbs, ubs = bounds(mm)
    ocs = objective(mm)

    json = Dict{String,Any}()
    json["id"] = "model" # default

    json[first(constants.keynames.genes)] = [
        Dict([
            "id" => gid,
            "name" => gene_name(mm, gid),
            "annotation" => gene_annotations(mm, gid),
            "notes" => gene_notes(mm, gid),
        ],) for gid in gene_ids
    ]

    json[first(constants.keynames.mets)] = [
        Dict([
            "id" => mid,
            "name" => metabolite_name(mm, mid),
            "formula" => maybemap(unparse_formula, metabolite_formula(mm, mid)),
            "charge" => metabolite_charge(mm, mid),
            "compartment" => metabolite_compartment(mm, mid),
            "annotation" => metabolite_annotations(mm, mid),
            "notes" => metabolite_notes(mm, mid),
        ]) for mid in met_ids
    ]

    json[first(constants.keynames.rxns)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            res["name"] = reaction_name(mm, rid)
            res["subsystem"] = reaction_subsystem(mm, rid)
            res["annotation"] = reaction_annotations(mm, rid)
            res["notes"] = reaction_notes(mm, rid)

            grr = reaction_gene_association(mm, rid)
            if !isnothing(grr)
                res["gene_reaction_rule"] = unparse_grr(String, grr)
            end

            res["lower_bound"] = lbs[ri]
            res["upper_bound"] = ubs[ri]
            res["objective_coefficient"] = ocs[ri]
            I, V = findnz(S[:, ri])
            res["metabolites"] =
                Dict{String,Float64}([met_ids[ii] => vv for (ii, vv) in zip(I, V)])
            res
        end for (ri, rid) in enumerate(rxn_ids)
    ]

    return JSONModel(json)
end
