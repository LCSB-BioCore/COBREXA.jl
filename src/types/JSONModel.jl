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

### Generic interface

# Unfortunately, this model type does not have standardized field names, hence the need to look for valid field names.
# The keys used to look for valid field names are in `src/base/constants.jl`.

"""
    _get_reactions(model::JSONModel)

Return a list of reaction dicts from `JSONModel`.
"""
function _get_reactions(model::JSONModel)
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            @info "Used key: \"$k\" to access reactions."
            return model.m[k][:]
        end
    end
    @warn "No reactions found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

"""
    _get_metabolites(model::JSONModel)

Return a list of metabolite dicts from `JSONModel`.
"""
function _get_metabolites(model::JSONModel)
    for k in _constants.possible_met_keys
        if haskey(model.m, k)
            @info "Used key: \"$k\" to access metabolites."
            return model.m[k][:]
        end
    end
    @warn "No metabolites found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

"""
    _get_genes(model::JSONModel)

Return a list of gene dicts from `JSONModel`.
"""
function _get_genes(model::JSONModel)
    for k in _constants.possible_gene_keys
        if haskey(model.m, k)
            @info "Used key: \"$k\" to access genes."
            return model.m[k][:]
        end
    end
    @warn "No genes found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function reactions(model::JSONModel)    
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"id\" is the key for reaction dicts..."
        return [string(get(r, "id", "")) for r in rxns] # assume `id` is the dict key.
    end
end

n_reactions(model::JSONModel) = length(reactions(model))

function metabolites(model::JSONModel)
    mets = _get_metabolites(model)
    if !isnothing(mets)
        @info "Assuming \"id\" is the key for metabolite dicts..."
        return [string(get(m, "id", "")) for m in mets] # assume `id` is the dict key.
    end
end

n_metabolites(model::JSONModel) = length(metabolites(model))

function genes(model::JSONModel)
    gs = _get_genes(model)
    if !isnothing(gs)
        @info "Assuming \"id\" is the key for gene dicts..."
        return [string(get(g, "id", "")) for g in gs] # assume `id` is the dict key.
    end
end

n_genes(model::JSONModel) = length(genes(model))

function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    rxn_key = _constants.possible_rxn_keys[[
        haskey(model.m, x) for x in _constants.possible_rxn_keys
    ]][1] # get the rxn key used
    met_ids = metabolites(model)
    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    @info "Reconstructing stoichiometric matrix..."
    for (i, rxn_id) in enumerate(rxn_ids)
        rxn_dict = model.m[rxn_key][i]["metabolites"] # assume metabolites is the only possible key
        for (met_id, coeff) in rxn_dict # assume met_id => coeff dict
            j = findfirst(x -> x == met_id, met_ids)
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
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"lower_bound\" is the key for reaction dicts... Otherwise default lower bound is used."
        for rxn in rxns
            push!(lbs, get(rxn, "lower_bound", -_constants.default_reaction_bound)) # assume `lower_bound` is the dict key
        end
        return lbs
    end
end

function upper_bounds(model::JSONModel)
    ubs = Float64[]
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"upper_bound\" is the key for reaction dicts... Otherwise default upper bound is used."
        for rxn in rxns
            push!(ubs, get(rxn, "upper_bound", _constants.default_reaction_bound)) # assume `upper_bound` is the dict key
        end
        return ubs
    end
end

function bounds(model::JSONModel)
    return lower_bounds(model), upper_bounds(model)
end

function balance(model::JSONModel)
    @info "Assuming balance is identically zero..."
    return spzeros(n_metabolites(model))
end

function objective(model::JSONModel)
    c = spzeros(n_reactions(model))
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"objective_coefficient\" is the key used for objective reactions..."
        for (i, rxn) in enumerate(rxns)
            if haskey(rxn, "objective_coefficient") # assume the only key?
                c[i] = rxn["objective_coefficient"]
            end
        end
        nnz(c) == 0 && (@warn "No objective found.")
        return c
    end
end

function gene_associations(model::JSONModel)
    grrs = String[]
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"gene_reaction_rule\" is the key used for gene reaction rules in reactions..."
        for rxn in rxns
            push!(grrs, get(rxn, "gene_reaction_rule", "")) # assume only key
        end
        length(grrs) == 0 && (@warn "No GRRs found.")
        return grrs
    end
end

function formulas(model::JSONModel)
    formula_strings = String[]
    mets = _get_metabolites(model)
    if !isnothing(mets)
        @info "Assuming \"formula\" is the key used for formulas in metabolites..."
        for met in mets
            push!(formula_strings, get(met, "formula", "")) # assume only key
        end
    end
    return formula_strings
end

function charges(model::JSONModel)
    charges_arr = Int64[] # assume only integer charges :) 
    mets = _get_metabolites(model)
    if !isnothing(mets)
        @info "Assuming \"charge\" is the key used for charges in metabolites..."
        for met in mets
            push!(charges_arr, get(met, "charge", 0)) # assume only key
        end
    end
    return charges_arr
end

function metabolite_chemistry(model::JSONModel)
    fs = formulas(model)
    cs = charges(model)
    return fs, cs
end

function reaction_subsystems(model::JSONModel)
    rsub = String[]
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        @info "Assuming \"subsystem\" is the key used for subsystems in reactions..."
        for rxn in rxns
            push!(rsub, get(rxn, "subsystem", "")) # assume only key
        end
        return rsub
    end
end

function metabolite_compartments(model::JSONModel)
    compartments = String[]
    mets = _get_metabolites(model)
    if !isnothing(mets)
        @info "Assuming \"compartment\" is the key used for compartment in metabolites..."
        for met in mets
            push!(compartments, get(met, "compartment", "")) # assume only key
        end
    end
    return compartments
end

function metabolite_notes(model::JSONModel)
    mets = _get_metabolites(model)
    if !isnothing(mets)
        m_notes = Vector{Any}() #Vector{Dict{String, Vector{String}}}()
        for m in mets
            push!(m_notes, _notes_from_jsonmodel(m))
        end
        return m_notes
    end
    return nothing
end

function metabolite_annotations(model::JSONModel)
    mets = _get_metabolites(model)
    if !isnothing(mets)
        m_annos = Vector{Dict{String, Vector{String}}}()
        for m in mets
            push!(m_annos, _annotation_from_jsonmodel(m))
        end
        return m_annos 
    end
    return nothing
end

function gene_notes(model::JSONModel)
    gs = _get_genes(model)
    if !isnothing(gs)
        g_notes = Vector{Dict{String, Vector{String}}}()
        for g in gs
            push!(g_notes, _notes_from_jsonmodel(g))
        end
        return g_notes
    end
    return nothing
end

function gene_annotations(model::JSONModel)
    gs = _get_genes(model)
    if !isnothing(gs)
        g_annos = Vector{Dict{String, Vector{String}}}()
        for g in gs
            push!(g_annos, _annotation_from_jsonmodel(g))
        end
        return g_annos
    end
    return nothing
end

function reaction_notes(model::JSONModel)
    rxns = _get_metabolites(model)
    if !isnothing(rxns)
        r_notes = Vector{Dict{String, Vector{String}}}()
        for r in rxns
            push!(r_notes, _notes_from_jsonmodel(r))
        end
        return r_notes
    end
    return nothing
end

function reaction_annotations(model::JSONModel)
    rxns = _get_metabolites(model)
    if !isnothing(rxns)
        r_annos = Vector{Dict{String, Vector{String}}}()
        for r in rxns
            push!(r_annos, _annotation_from_jsonmodel(r))
        end
        return r_annos
    end
    return nothing
end

### Accessor functions to construct StandardModel from a JSONModel efficiently (loop through model struct only once)

function _gene_ordereddict(model::JSONModel)
    gs = _get_genes(model)
    if !isnothing(gs)
        gd = OrderedDict{String,Gene}()
        for g in gs
            gg = Gene(
                g["id"];
                name = get(g, "name", ""),
                notes = _notes_from_jsonmodel(g),
                annotation = _annotation_from_jsonmodel(g),
            )
            gd[g["id"]] = gg
        end
        return gd
    end
end

function _reaction_ordereddict(model::JSONModel)
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        rd = OrderedDict{String,Reaction}()
        for r in rxns
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
                notes = _notes_from_jsonmodel(r),
                annotation = _annotation_from_jsonmodel(r),
                objective_coefficient = get(r, "objective_coefficient", 0.0),
            )
            rd[r["id"]] = rr
        end
        return rd
    end
end

function _metabolite_ordereddict(model::JSONModel)
    mets = _get_metabolites(model)
    if !isnothing(mets)
        md = OrderedDict{String,Metabolite}()
        for m in mets
            mm = Metabolite(
                m["id"];
                name = get(m, "name", ""),
                formula = get(m, "formula", ""),
                charge = get(m, "charge", 0),
                compartment = get(m, "compartment", ""),
                notes = _notes_from_jsonmodel(m),
                annotation = _annotation_from_jsonmodel(m),
            )
            md[m["id"]] = mm
        end
       return md
    end
end

# convert to Dict{String, Vector{String}}
function _notes_from_jsonmodel(x; verbose=false)
    verbose && @info "Assuming \"notes\" is the key used for storing notes"
    d = get(x, "notes", Dict{String,Vector{String}}())
    dd = Dict{String,Vector{String}}()
    for (k, v) in d
        dd[k] = string.(v)
    end
    return dd
end

# convert to Dict{String, Vector{String}}
function _annotation_from_jsonmodel(x; verbose=false)
    verbose && @info "Assuming \"annotation\" is the key used for storing annotations"
    d = get(x, "annotation", Dict{String,Union{Vector{String},String}}())
    dd = Dict{String, Vector{String}}()
    for (k, v) in d
        if k == "sbo" || k == "SBO" # sbo terms are not assigned to arrays in JSON models
            dd[k] = [string(v)]
        else
            dd[k] = string.(v)
        end
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
