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
