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

    return sparse([get(rxn, "objective_coefficient", 0) for rxn in model.m[r]])
end

"""
    reaction_gene_associaton(model::JSONModel, rid::String)

Parses the `.gene_reaction_rule` from reactions.
"""
function reaction_gene_associaton(model::JSONModel, rid::String)
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

#
# Below lies the batch-processing getter API, arguably more efficient for
# making StdModel from JSON directly.
#

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
        m_annos = Vector{Dict{String,Vector{String}}}()
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
        g_notes = Vector{Dict{String,Vector{String}}}()
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
        g_annos = Vector{Dict{String,Vector{String}}}()
        for g in gs
            push!(g_annos, _annotation_from_jsonmodel(g))
        end
        return g_annos
    end
    return nothing
end

function reaction_notes(model::JSONModel)
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        r_notes = Vector{Dict{String,Vector{String}}}()
        for r in rxns
            push!(r_notes, _notes_from_jsonmodel(r))
        end
        return r_notes
    end
    return nothing
end

function reaction_annotations(model::JSONModel)
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        r_annos = Vector{Dict{String,Vector{String}}}()
        for r in rxns
            push!(r_annos, _annotation_from_jsonmodel(r))
        end
        return r_annos
    end
    return nothing
end

function reaction_names(model::JSONModel)
    rxns = _get_reactions(model)
    if !isnothing(rxns)
        r_names = String[]
        @info "Assuming \"name\" is the key used to get the reaction names..."
        for r in rxns
            push!(r_names, get(r, "name", ""))
        end
        return r_names
    end
    return nothing
end

function metabolite_names(model::JSONModel)
    mets = _get_metabolites(model)
    if !isnothing(mets)
        m_names = String[]
        @info "Assuming \"name\" is the key used to get the metabolite names..."
        for m in mets
            push!(m_names, get(m, "name", ""))
        end
        return m_names
    end
    return nothing
end

function gene_names(model::JSONModel)
    gs = _get_genes(model)
    if !isnothing(gs)
        g_names = String[]
        @info "Assuming \"name\" is the key used to get the gene names..."
        for g in gs
            push!(g_names, get(g, "name", ""))
        end
        return g_names
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
function _notes_from_jsonmodel(x; verbose = false)
    verbose && @info "Assuming \"notes\" is the key used for storing notes"
    d = get(x, "notes", Dict{String,Vector{String}}())
    dd = Dict{String,Vector{String}}()
    for (k, v) in d
        dd[k] = string.(v)
    end
    return dd
end

# convert to Dict{String, Vector{String}}
function _annotation_from_jsonmodel(x; verbose = false)
    verbose && @info "Assuming \"annotation\" is the key used for storing annotations"
    d = get(x, "annotation", Dict{String,Union{Vector{String},String}}())
    dd = Dict{String,Vector{String}}()
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

function Base.convert(::typeof{JSONModel}, mm::MetabolicModel)
    rxn_ids = reactions(mm)
    met_ids =  metabolites(mm)
    gene_ids = genes(mm)
    S = stoichiometry(mm)
    lbs, ubs = bounds(mm)
    grrs = gene_annotations(mm)
    
    r_annos = reaction_annotations(mm)
    r_notes = reaction_notes(mm)
    m_annos = metabolite_annotations(mm)
    m_notes = metabolite_notes(mm)
    g_annos = gene_annotations(mm)
    g_notes = gene_notes(mm)
    
    fs, cs = metabolite_chemistry(mm) # formulas, charges
    r_subs = reaction_subsystems(mm)
    m_comps = metabolite_compartments(mm)

    objs = objective(mm)

    r_names = reaction_names(model)
    g_names = gene_names(model)
    m_names = metabolite_names(model)

    json_model = JSONModel(Dict{String,Any}())
    json_model.m["id"] = id(mm)

    # build genes
    gs = Any[]
    for i=1:n_genes(mm)
        g = Dict{String, Any}()
        g["id"] = gene_ids[i]
        !isnothing(g_names) && (g["name"] = g_names[i])
        !isnothing(g_notes) && (g["notes"] = g_notes[i])
        !isnothing(g_annos) && (g["annotation"] = g_annos[i])
        push!(gs, g)
    end
    json_model.m["genes"] = gs

    # build metabolites
    ms = Any[]
    for i=1:n_metabolites(mm)
        m = Dict{String, Any}()
        m["id"] = met_ids[i]
        !isnothing(m_comps) && (m["compartment"] = m_comps[i])
        !isnothing(fs) && (m["formula"] = fs[i])
        !isnothing(cs) && (m["charge"] = cs[i])
        !isnothing(m_notes) && (m["notes"] = m_notes[i])
        !isnothing(m_annos) && (m["annotation"] = m_annos[i])
        push!(ms, m)
    end
    json_model.m["metabolites"] = ms

    # build reactions
    rs = Any[]
    for i=1:n_reactions(mm)
        r = Dict{String, Any}()
        r["id"] = rxn_ids[i]
        !isnothing(r_notes) && (r["notes"] = r_notes[i])
        !isnothing(r_annos) && (r["annotation"] = r_annos[i])
        !isnothing(r_subs) && (r["subsystem"] = r_subs[i])
        !isnothing(grrs) && (r["gene_reaction_rule"] = grrs[i])
        !isnothing(ubs) && (r["upper_bound"] = ubs[i])
        !isnothing(lbs) && (r["lower_bound"] = lbs[i])
        !isnothing(r_names) && (r["name"] = r_names[i])
        
        reqn = Dict{String, Float64}()
        I, V = findnz(S[:, i])
        for (ii, vv) in zip(I, V)
            reqn[met_ids[ii]] = vv
        end
        r["metabolites"] = reqn
        
        push!(rs, r)
    end
    json_model.m["reactions"] = rs

    return json_model
end
