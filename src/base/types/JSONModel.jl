"""
    struct JSONModel

A struct used to store the contents of a JSON model, i.e. a model read from a file ending with `.json`. 
These model files typically store all the model parameters in arrays of dictionaries. 

When importing to this model type, no information gets lost between the file and the object in memory.
However, not all of the fields can be used in analysis functions, and not all fields are captured when converting to `StandardModel`

Note, this model type is not very efficient, especially when calling many generic interface functions sequentially.
Instead use one of the COBREXA model types: `StandardModel`, `CoreModel` or `CoreModelCoupled` if speed is important.

See also: [`CoreModel`](@ref), [`CoreModelCoupled`](@ref), [`StandardModel`](@ref)

# Example
````
model = load_json_model("some_model.json")
model.m # the actual underlying JSON
````
"""
struct JSONModel <: MetabolicModel
    m::Dict{String,Any}
end


# helper access macros, see examples below for usage
macro _json_sectionkey(namesConstant::Symbol, error)
    esc(:(
        begin
            _key = _guesskey(keys(model.m), _constants.keynames.$namesConstant)
            if isnothing(_key)
                $error
            end
            _key
        end
    ))
end

macro _json_section(namesConstant::Symbol, error)
    esc(:(
        begin
            _key = @_json_sectionkey $namesConstant ($error)
            model.m[_key]
        end
    ))
end

macro _json_section_firstid(namesConstant::Symbol, id::Symbol, error)
    return esc(:(
        begin
            _s = @_json_section $namesConstant ($error)
            _firstmatch(x -> x["id"] == $id, _s)
        end
    ))
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
function reactions(model::JSONModel)
    rs = @_json_section rxns throw(
        DomainError(keys(model.m), "JSON model has no reaction keys"),
    )

    return [string(get(rs[i], "id", "rxn$i")) for i in eachindex(rs)]
end

"""
    metabolites(model::JSONModel)

Extract metabolite names (stored as `.id`) from JSON model.
"""
function metabolites(model::JSONModel)
    ms = @_json_section mets throw(
        DomainError(keys(model.m), "JSON model has no metabolite keys"),
    )

    return [string(get(ms[i], "id", "met$i")) for i in eachindex(ms)]
end

"""
    genes(model::JSONModel)

Extract gene names from a JSON model.
"""
function genes(model::JSONModel)
    gs = @_json_section genes return []

    return [string(get(gs[i], "id", "gene$i")) for i in eachindex(gs)]
end

"""
    stoichiometry(model::JSONModel)

Get the stoichiometry. Assuming the information is stored in reaction object
under key `.metabolites`.
"""
function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    met_ids = metabolites(model)

    rs = @_json_section rxns throw(
        DomainError(keys(model.m), "JSON model has no reaction keys"),
    )

    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    for (i, rid) in enumerate(rxn_ids)
        r = _firstmatch(x -> x["id"] == rid, rs)
        for (met_id, coeff) in r["metabolites"]
            j = findfirst(==(met_id), met_ids)
            if isnothing(j)
                throw(
                    DomainError(
                        met_id,
                        "Unknown metabolite found in stoichiometry of $(rxn_ids[i])",
                    ),
                )
            end
            S[j, i] = coeff
        end
    end
    return S
end

"""
    bounds(model::JSONModel)

Get the bounds for reactions, assuming the information is stored in
`.lower_bound` and `.upper_bound`.
"""
function bounds(model::JSONModel)
    rs = @_json_section rxns return (
        sparse(fill(-_constants.default_reaction_bound, n_reactions(model))),
        sparse(fill(_constants.default_reaction_bound, n_reactions(model))),
    )
    return (
        sparse([get(rxn, "lower_bound", -_constants.default_reaction_bound) for rxn in rs]),
        sparse([get(rxn, "upper_bound", _constants.default_reaction_bound) for rxn in rs]),
    )
end

"""
    objective(model::JSONModel)

Collect `.objective_coefficient` keys from model reactions.
"""
function objective(model::JSONModel)
    rs = @_json_section rxns return spzeros(n_reactions(model))

    return sparse([float(get(rxn, "objective_coefficient", 0.0)) for rxn in rs])
end

"""
    reaction_gene_associaton(model::JSONModel, rid::String)

Parses the `.gene_reaction_rule` from reactions.
"""
function reaction_gene_association(model::JSONModel, rid::String)
    rxn = @_json_section_firstid rxns rid return nothing
    return maybemap(_parse_grr, get(rxn, "gene_reaction_rule", nothing))
end

"""
    reaction_subsystem(model::JSONModel, rid::String)

Parses the `.subsystem` out from reactions.
"""
function reaction_subsystem(model::JSONModel, rid::String)
    rxn = @_json_section_firstid rxns rid return nothing
    return get(rxn, "subsystem", nothing)
end

"""
    metabolite_formula(model::JSONModel, mid::String)

Parse and return the metabolite `.formula`
"""
function metabolite_formula(model::JSONModel, mid::String)
    met = @_json_section_firstid mets mid return nothing
    return maybemap(_formula_to_atoms, get(met, "formula", nothing))
end

"""
    metabolite_charge(model::JSONModel, mid::String)

Return the metabolite `.charge`
"""
function metabolite_charge(model::JSONModel, mid::String)
    met = @_json_section_firstid mets mid return nothing
    return get(met, "charge", 0)
end

"""
    metabolite_compartment(model::JSONModel, mid::String)

Return the metabolite `.compartment`
"""
function metabolite_compartment(model::JSONModel, mid::String)
    met = @_json_section_firstid mets mid return nothing
    return get(met, "compartment", nothing)
end

function gene_annotations(model::JSONModel, gid::String)::Annotations
    gene = @_json_section_firstid genes gid return Dict()
    maybemap(_parse_annotations, get(gene, "annotation", nothing))
end

function gene_notes(model::JSONModel, gid::String)::Notes
    gene = @_json_section_firstid genes gid return Dict()
    maybemap(_parse_notes, get(gene, "notes", nothing))
end

function reaction_annotations(model::JSONModel, rid::String)::Annotations
    rxn = @_json_section_firstid rxns rid return Dict()
    maybemap(_parse_annotations, get(rxn, "annotation", nothing))
end

function reaction_notes(model::JSONModel, rid::String)::Notes
    rxn = @_json_section_firstid rxns rid return Dict()
    maybemap(_parse_notes, get(rxn, "notes", nothing))
end

function metabolite_annotations(model::JSONModel, mid::String)::Annotations
    met = @_json_section_firstid mets mid return Dict()
    maybemap(_parse_annotations, get(met, "annotation", nothing))
end

function metabolite_notes(model::JSONModel, mid::String)::Notes
    met = @_json_section_firstid mets mid return Dict()
    maybemap(_parse_notes, get(met, "notes", nothing))
end


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

    json_model = JSONModel(Dict{String,Any}())
    m = Dict{String,Any}()
    m["id"] = "model" # default

    #TODO: add notes, names and similar fun stuff when they are available

    m[first(_constants.keynames.genes)] = [
        Dict([
            "id" => gid,
            "annotation" => gene_annotations(mm, gid),
            "notes" => gene_notes(mm, gid),
        ],) for gid in gene_ids
    ]

    m[first(_constants.keynames.mets)] = [
        Dict([
            "id" => mid,
            "formula" => maybemap(_atoms_to_formula, metabolite_formula(mm, mid)),
            "charge" => metabolite_charge(mm, mid),
            "compartment" => metabolite_compartment(mm, mid),
            "annotation" => metabolite_annotations(mm, mid),
            "notes" => metabolite_notes(mm, mid),
        ]) for mid in met_ids
    ]

    m[first(_constants.keynames.rxns)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            res["subsystem"] = reaction_subsystem(mm, rid)
            res["annotation"] = reaction_annotations(mm, rid)
            res["notes"] = reaction_notes(mm, rid)

            grr = reaction_gene_association(mm, rid)
            if !isnothing(grr)
                res["gene_reaction_rule"] = _unparse_grr(grr)
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

    return JSONModel(m)
end
