"""
    struct JSONModel <: MetabolicModel
        json::Dict{String,Any}
    end

A struct used to store the contents of a JSON model, i.e. a model read from a
file ending with `.json`. These model files typically store all the model
parameters in arrays of JSON objects (i.e. Julia dictionaries).

Usually, not all of the fields of the input JSON can be easily represented when
converting to other models, care should be taken to avoid losing information.

Direct work on this precise model type is not very efficient, as the accessor
functions need to repeatedly find the information in the JSON tree. This gets
very slow especially if calling many accessor functions sequentially. To avoid
that, convert to e.g. [`StandardModel`](@ref) as soon as possible.

# Example
````
model = load_json_model("some_model.json")
model.json # see the actual underlying JSON
reactions(model) # see the list of reactions
````
"""
struct JSONModel <: MetabolicModel
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
    rkey = _guesskey(keys(json), _constants.keynames.rxns)
    isnothing(rkey) && throw(DomainError(keys(json), "JSON model has no reaction keys"))
    rs = json[rkey]

    mkey = _guesskey(keys(json), _constants.keynames.mets)
    ms = json[mkey]
    isnothing(mkey) && throw(DomainError(keys(json), "JSON model has no metabolite keys"))

    gkey = _guesskey(keys(json), _constants.keynames.genes)
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

"""
    reactions(model::JSONModel)

Extract reaction names (stored as `.id`) from JSON model.
"""
reactions(model::JSONModel) = [_json_rxn_name(r, i) for (i, r) in enumerate(model.rxns)]

"""
    metabolites(model::JSONModel)

Extract metabolite names (stored as `.id`) from JSON model.
"""
metabolites(model::JSONModel) = [_json_met_name(m, i) for (i, m) in enumerate(model.mets)]

"""
    genes(model::JSONModel)

Extract gene names from a JSON model.
"""
genes(model::JSONModel) = [_json_gene_name(g, i) for (i, g) in enumerate(model.genes)]

"""
    stoichiometry(model::JSONModel)

Get the stoichiometry. Assuming the information is stored in reaction object
under key `.metabolites`.
"""
function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    met_ids = metabolites(model)

    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    for (i, rid) in enumerate(rxn_ids)
        r = model.rxns[model.rxn_index[rid]]
        for (mid, coeff) in r["metabolites"]
            if !haskey(model.met_index, mid)
                throw(
                    DomainError(
                        met_id,
                        "Unknown metabolite found in stoichiometry of $(rxn_ids[i])",
                    ),
                )
            end
            S[model.met_index[mid], i] = coeff
        end
    end
    return S
end

"""
    bounds(model::JSONModel)

Get the bounds for reactions, assuming the information is stored in
`.lower_bound` and `.upper_bound`.
"""
bounds(model::JSONModel) = (
    sparse([
        get(rxn, "lower_bound", -_constants.default_reaction_bound) for rxn in model.rxns
    ]),
    sparse([
        get(rxn, "upper_bound", _constants.default_reaction_bound) for rxn in model.rxns
    ]),
)

"""
    objective(model::JSONModel)

Collect `.objective_coefficient` keys from model reactions.
"""
objective(model::JSONModel) =
    sparse([float(get(rxn, "objective_coefficient", 0.0)) for rxn in model.rxns])

"""
    reaction_gene_associaton(model::JSONModel, rid::String)

Parses the `.gene_reaction_rule` from reactions.
"""
reaction_gene_association(model::JSONModel, rid::String) = _maybemap(
    _parse_grr,
    get(model.rxns[model.rxn_index[rid]], "gene_reaction_rule", nothing),
)

"""
    reaction_subsystem(model::JSONModel, rid::String)

Parses the `.subsystem` out from reactions.
"""
reaction_subsystem(model::JSONModel, rid::String) =
    get(model.rxns[model.rxn_index[rid]], "subsystem", nothing)

"""
    metabolite_formula(model::JSONModel, mid::String)

Parse and return the metabolite `.formula`
"""
metabolite_formula(model::JSONModel, mid::String) =
    _maybemap(_parse_formula, get(model.mets[model.met_index[mid]], "formula", nothing))

"""
    metabolite_charge(model::JSONModel, mid::String)

Return the metabolite `.charge`
"""
metabolite_charge(model::JSONModel, mid::String) =
    get(model.mets[model.met_index[mid]], "charge", 0)

"""
    metabolite_compartment(model::JSONModel, mid::String)

Return the metabolite `.compartment`
"""
metabolite_compartment(model::JSONModel, mid::String) =
    get(model.mets[model.met_index[mid]], "compartment", nothing)

"""
    gene_annotations(model::JSONModel, gid::String)::Annotations

Gene annotations from the [`JSONModel`](@ref).
"""
gene_annotations(model::JSONModel, gid::String)::Annotations = _maybemap(
    _parse_annotations,
    get(model.genes[model.gene_index[gid]], "annotation", nothing),
)

"""
    gene_notes(model::JSONModel, gid::String)::Notes

Gene notes from the [`JSONModel`](@ref).
"""
gene_notes(model::JSONModel, gid::String)::Notes =
    _maybemap(_parse_notes, get(model.genes[model.gene_index[gid]], "notes", nothing))

"""
    reaction_annotations(model::JSONModel, rid::String)::Annotations

Reaction annotations from the [`JSONModel`](@ref).
"""
reaction_annotations(model::JSONModel, rid::String)::Annotations = _maybemap(
    _parse_annotations,
    get(model.rxns[model.rxn_index[rid]], "annotation", nothing),
)

"""
    reaction_notes(model::JSONModel, rid::String)::Notes

Reaction notes from the [`JSONModel`](@ref).
"""
reaction_notes(model::JSONModel, rid::String)::Notes =
    _maybemap(_parse_notes, get(model.rxns[model.rxn_index[rid]], "notes", nothing))

"""
    metabolite_annotations(model::JSONModel, mid::String)::Annotations

Metabolite annotations from the [`JSONModel`](@ref).
"""
metabolite_annotations(model::JSONModel, mid::String)::Annotations = _maybemap(
    _parse_annotations,
    get(model.mets[model.met_index[mid]], "annotation", nothing),
)

"""
    metabolite_notes(model::JSONModel, mid::String)::Notes

Metabolite notes from the [`JSONModel`](@ref).
"""
metabolite_notes(model::JSONModel, mid::String)::Notes =
    _maybemap(_parse_notes, get(model.mets[model.met_index[mid]], "notes", nothing))

"""
    reaction_stoichiometry(model::JSONModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(model::JSONModel, rid::String)::Dict{String,Float64} =
    model.rxns[model.rxn_index[rid]]["metabolites"]

"""
    Base.convert(::Type{JSONModel}, mm::MetabolicModel)

Convert any [`MetabolicModel`](@ref) to [`JSONModel`](@ref).
"""
function Base.convert(::Type{JSONModel}, mm::MetabolicModel)
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

    #TODO: add notes, names and similar fun stuff when they are available

    json[first(_constants.keynames.genes)] = [
        Dict([
            "id" => gid,
            "annotation" => gene_annotations(mm, gid),
            "notes" => gene_notes(mm, gid),
        ],) for gid in gene_ids
    ]

    json[first(_constants.keynames.mets)] = [
        Dict([
            "id" => mid,
            "formula" => _maybemap(_unparse_formula, metabolite_formula(mm, mid)),
            "charge" => metabolite_charge(mm, mid),
            "compartment" => metabolite_compartment(mm, mid),
            "annotation" => metabolite_annotations(mm, mid),
            "notes" => metabolite_notes(mm, mid),
        ]) for mid in met_ids
    ]

    json[first(_constants.keynames.rxns)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            res["subsystem"] = reaction_subsystem(mm, rid)
            res["annotation"] = reaction_annotations(mm, rid)
            res["notes"] = reaction_notes(mm, rid)

            grr = reaction_gene_association(mm, rid)
            if !isnothing(grr)
                res["gene_reaction_rule"] = _unparse_grr(String, grr)
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
