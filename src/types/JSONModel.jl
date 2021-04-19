struct JSONModel <: MetabolicModel
    m::Dict{String,Any}
end

### Generic interface
# Unfortunately this model type does not have standardized field names, hence the need to look for valid fieldnames.
# The keys used to look for valid fieldnames is in `constants`.
# However, assume that `id` is used to index reactions, metabolites and genes.
# Also assume that `metabolites` is used to access the dict containing the reaction equation for a reaction.

function reactions(model::JSONModel)
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            return [string(r["id"]) for r in model.m[k][:]]
        end
    end
    @warn "No reactions found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_reactions(model::JSONModel) = length(reactions(model))

function metabolites(model::JSONModel)
    for k in _constants.possible_met_keys
        if haskey(model.m, k)
            return [string(m["id"]) for m in model.m[k][:]]
        end
    end
    @warn "No metabolites found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_metabolites(model::JSONModel) = length(metabolites(model))

function genes(model::JSONModel)
    for k in _constants.possible_gene_keys
        if haskey(model.m, k)
            return [string(g["id"]) for g in model.m[k][:]]
        end
    end
    @warn "No genes found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_genes(model::JSONModel) = length(genes(model))

function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    rxn_key = _constants.possible_rxn_keys[[
        haskey(model.m, x) for x in _constants.possible_rxn_keys
    ]][1] # get the rxn key used
    met_ids = metabolites(model)
    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    for (i, rxn_id) in enumerate(rxn_ids)
        rxn_dict = model.m[rxn_key][i]["metabolites"] # assume metabolites is the only possible key
        for (met_id, coeff) in rxn_dict # assume met_id => coeff dict
            j = findfirst(x -> x == met_id, met_ids) # row
            isnothing(j) ?
            (@error "S matrix construction error: $(met_id) not defined."; return nothing) :
            nothing
            S[j, i] = coeff
        end
    end
    return S
end

function lower_bounds(model::JSONModel)
    lbs = Float64[]
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for rxn in model.m[k]
                push!(lbs, rxn["lower_bound"]) # assume this is the only possible key
            end
            return lbs
        end
    end
    @warn "No lower bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function upper_bounds(model::JSONModel)
    ubs = Float64[]
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for rxn in model.m[k]
                push!(ubs, rxn["upper_bound"]) # assume this is the only possible key
            end
            return ubs
        end
    end
    @warn "No upper bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function bounds(model::JSONModel)
    return lower_bounds(model), upper_bounds(model)
end

function balance(model::JSONModel)
    return spzeros(n_metabolites(model))
end

function objective(model::JSONModel)
    c = spzeros(n_reactions(model))
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for (i, rxn) in enumerate(model.m[k])
                if haskey(rxn, "objective_coefficient") # assume the only key?
                    c[i] = rxn["objective_coefficient"]
                end
            end
            return c
        end
    end
    return nothing
end

function gene_reaction_rules(model::JSONModel)
    grrs = String[]
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for rxn in model.m[k]
                push!(grrs, rxn["gene_reaction_rules"]) # assume only key
            end
            return grrs
        end
    end
    return nothing
end

function formulas(model::JSONModel)
    formula_strings = String[]
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for rxn in model.m[k]
                push!(grrs, rxn["gene_reaction_rules"]) # assume only key
            end
        end
    end
    return formula_strings
end

### Accessor functions to construct StandardModel efficiently (loop through model struct only once)

function _gene_ordereddict(model::JSONModel)
    gd = OrderedDict{String,Gene}()
    for k in _constants.possible_gene_keys
        if haskey(model.m, k)
            for g in model.m[k]
                gg = Gene(
                    g["id"];
                    name = get(g, "name", ""),
                    notes = _notes_from_jsonmodel(
                        get(g, "notes", Dict{String,Vector{String}}()),
                    ),
                    annotation = _annotation_from_jsonmodel(
                        get(g, "annotation", Dict{String,Union{Vector{String},String}}()),
                    ),
                )
                gd[g["id"]] = gg
            end
            break
        end
    end
    return gd
end

function _reaction_ordereddict(model::JSONModel)
    rd = OrderedDict{String,Reaction}()
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for r in model.m[k]
                rr = Reaction(
                    r["id"];
                    name = get(r, "name", ""),
                    metabolites = _reaction_formula_from_jsonmodel(
                        get(r, "metabolites", Dict{String,Any}()),
                    ),
                    lb = get(r, "lb", -_constants.default_reaction_bound),
                    ub = get(r, "ub", _constants.default_reaction_bound),
                    grr = _parse_grr(get(r, "gene_reaction_rule", "")),
                    subsystem = get(r, "subsystem", ""),
                    notes = _notes_from_jsonmodel(get(r, "notes", Dict{String,Any}())),
                    annotation = _annotation_from_jsonmodel(
                        get(r, "annotation", Dict{String,Union{Vector{String},String}}()),
                    ),
                    objective_coefficient = get(r, "objective_coefficient", 0.0),
                )
                rd[r["id"]] = rr
            end
            break
        end
    end
    return rd
end

function _metabolite_ordereddict(model::JSONModel)
    md = OrderedDict{String,Metabolite}()
    for k in _constants.possible_met_keys
        if haskey(model.m, k)
            for m in model.m[k]
                mm = Metabolite(
                    m["id"];
                    name = get(m, "name", ""),
                    formula = get(m, "formula", ""),
                    charge = get(m, "charge", 0),
                    compartment = get(m, "compartment", ""),
                    notes = _notes_from_jsonmodel(
                        get(m, "notes", Dict{String,Vector{String}}()),
                    ),
                    annotation = _annotation_from_jsonmodel(
                        get(m, "annotation", Dict{String,Union{Vector{String},String}}()),
                    ),
                )
                md[m["id"]] = mm
            end
            break
        end
    end
    return md
end

# convert d to Dict{String, Vector{String}}
function _notes_from_jsonmodel(d)
    dd = Dict{String,Vector{String}}()
    for (k, v) in d
        dd[k] = string.(v)
    end
    return dd
end

# convert d to Dict{String, Union{String, Vector{String}}}
function _annotation_from_jsonmodel(d)
    dd = Dict{String,Union{String,Vector{String}}}()
    for (k, v) in d
        dd[k] = string.(v)
    end
    return dd
end

# convert d to Dict{String, Float64}
function _reaction_formula_from_jsonmodel(d)
    dd = Dict{String,Float64}()
    for (k, v) in d
        dd[k] = float(v)
    end
    return dd
end

# Construct a JSONModel from any other model

function Base.convert(::typeof(JSONModel), model::MetabolicModel)
    jsonmodel = JSONModel(Dict{String,Any}()) # blank model
    met_ids = metabolites(model)
    rxn_ids = reactions(model)
    gene_ids = genes(model)
    grrs = gene_reaction_rules(model)
    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    objs = objective(model)

end
