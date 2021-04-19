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
model = read_model("some_model.json")
model.m # the actual underlying model
````
"""
struct JSONModel <: MetabolicModel
    m::Dict{String,Any}
end


function id(model::JSONModel)
    for k in _constants.keynames.ids
        if haskey(model.m, k)
            @info "Used key: \"$k\" to access the model id."
            return model.m[k]
        end
    end
end

"""
    _guesskey(ks, possibilities)

Unfortunately, JSON models do not have standardized field names, so we need to
try a few possibilities and guess the best one.  The keys used to look for
valid field names are specified in `src/base/constants.jl`.
"""
function _guesskey(avail, possibilities)
    x = intersect(possibilities, avail)

    if isempty(x)
        @debug "could not find any of keys: $possibilities"
        return nothing
    end

    if length(x) > 1
        @debug "Possible ambiguity between keys: $x"
    end
    return x[1]
end

"""
    reactions(model::JSONModel)

Extract reaction names (stored as `.id`) from JSON model.
"""
function reactions(model::JSONModel)
    k = _guesskey(keys(model.m), _constants.keynames.rxns)
    if isnothing(k)
        throw(DomainError(keys(model.m), "JSON model has no reaction keys"))
    end

    return [string(get(model.m[k][i], "id", "rxn$i")) for i in eachindex(model.m[k])]
end

"""
    metabolites(model::JSONModel)

Extract metabolite names (stored as `.id`) from JSON model.
"""
function metabolites(model::JSONModel)
    k = _guesskey(keys(model.m), _constants.keynames.mets)
    if isnothing(k)
        throw(DomainError(keys(model.m), "JSON model has no metabolite keys"))
    end

    return [string(get(model.m[k][i], "id", "met$i")) for i in eachindex(model.m[k])]
end

"""
    genes(model::JSONModel)

Extract gene names from a JSON model.
"""
function genes(model::JSONModel)
    k = _guesskey(keys(model.m), _constants.keynames.genes)
    if isnothing(k)
        return [] #no genes
    end

    return [string(get(model.m[k][i], "id", "gene$i")) for i in eachindex(model.m[k])]
end

"""
    stoichiometry(model::JSONModel)

Get the stoichiometry. Assuming the information is stored in reaction object
under key `.metabolites`.
"""
function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    met_ids = metabolites(model)

    r = _guesskey(keys(model.m), _constants.keynames.rxns)
    if isnothing(r)
        throw(DomainError(keys(model.m), "JSON model has no reaction keys"))
    end

    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    for i in eachindex(rxn_ids)
        for (met_id, coeff) in model.m[r][i]["metabolites"]
            j = findfirst(x -> x == met_id, met_ids)
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
    r = _guesskey(keys(model.m), _constants.keynames.rxns)
    if isnothing(r)
        return (
            sparse(fill(-_constants.default_reaction_bound, n_reactions(model))),
            sparse(fill(_constants.default_reaction_bound, n_reactions(model))),
        )
    end
    return (
        sparse([
            get(rxn, "lower_bound", -_constants.default_reaction_bound) for
            rxn in model.m[r]
        ]),
        sparse([
            get(rxn, "upper_bound", _constants.default_reaction_bound) for rxn in model.m[r]
        ]),
    )
end

"""
    objective(model::JSONModel)

Collect `.objective_coefficient` keys from model reactions.
"""
function objective(model::JSONModel)
    r = _guesskey(keys(model.m), _constants.keynames.rxns)
    if isnothing(r)
        return spzeros(n_reactions(model))
    end

    return sparse([float(get(rxn, "objective_coefficient", 0.0)) for rxn in model.m[r]])
end

"""
    reaction_gene_associaton(model::JSONModel, rid::String)

Parses the `.gene_reaction_rule` from reactions.
"""
function reaction_gene_association(model::JSONModel, rid::String)
    r = _guesskey(keys(model.m), _constants.keynames.rxns)
    if isnothing(r)
        return nothing
    end

    ri = first(indexin(rid, keys(model.m[r])))
    return maybemap(_parse_grr, get(model.m[r][ri], "gene_reaction_rule", nothing))
end

"""
    metabolite_chemistry(model::JSONModel, mid::String)

Parse and return the metabolite `.formula` and `.charge`.
"""
function metabolite_chemistry(model::JSONModel, mid::String)
    m = _guesskey(keys(model.m), _constants.keynames.mets)
    if isnothing(m)
        return nothing
    end

    mi = first(indexin(mid, keys(model.m[m])))
    met = models.m[m][mi]
    formula = maybemap(_formula_to_atoms, get(met, "formula", nothing))
    return maybemap(f -> (f, get(met, "charge", 0)), formula)
end

#TODO annotation accessors

function Base.convert(::Type{JSONModel}, mm::MetabolicModel)
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

    m[first(_constants.keynames.genes)] =
        [Dict(["id" => gid, "annotation" => gene_annotations(mm, gid)]) for gid in gene_ids]

    m[first(_constants.keynames.mets)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = mid
            ch = metabolite_chemistry(mm, mid)
            if !isnothing(ch)
                res["formula"] = _atoms_to_formula(ch[1])
                res["charge"] = ch[2]
            end
            res["annotation"] = metabolite_annotations(mm, mid)
            res
        end for mid in met_ids
    ]

    m[first(_constants.keynames.rxns)] = [
        begin
            res = Dict{String,Any}()
            res["id"] = rid
            a = reaction_annotations(mm, rid)
            if !isempty(a)
                res["annotation"] = a
            end
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
