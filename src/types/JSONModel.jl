struct JSONModel <: MetabolicModel
    m::Dict{String,Any}
end

# Generic interface
# Unfortunately this model type does not have standardized field names, hence the need to look for valid fieldnames.
# The keys used to look for valid fieldnames is in `constants`.
# However, assume that `id` is used to index reactions, metabolites and genes.
# Also assume that `metabolites` is used to access the dict containing the reaction equation for a reaction.

function reactions(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            return [string(r["id"]) for r in model.m[k][:]]
        end
    end
    @warn "No reactions found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_reactions(model::JSONModel)::Int = length(reactions(model))

function metabolites(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in _constants.possible_met_keys
        if haskey(model.m, k)
            return [string(m["id"]) for m in model.m[k][:]]
        end
    end
    @warn "No metabolites found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_metabolites(model::JSONModel)::Int = length(metabolites(model))

function genes(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in _constants.possible_gene_keys
        if haskey(model.m, k)
            return [string(g["id"]) for g in model.m[k][:]]
        end
    end
    @warn "No genes found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_genes(model::JSONModel)::Int = length(genes(model))

function stoichiometry(model::JSONModel)
    rxn_ids = reactions(model)
    rxn_key = _constants.possible_rxn_keys[[haskey(model.m, x) for x in _constants.possible_rxn_keys]][1] # get the rxn key used
    met_ids = metabolites(model)
    S = SparseArrays.spzeros(length(met_ids), length(rxn_ids))
    for (i, rxn_id) in enumerate(rxn_ids)
        rxn_dict = model.m[rxn_key][i]["metabolites"] # assume metabolites is the only possible key
        for (met_id, coeff) in rxn_dict # assume met_id => coeff dict
            j = findfirst(x -> x == met_id, met_ids) # row
            isnothing(j) ?
            (@error "S matrix construction error: $(met_id) not defined."; return nothing) : nothing
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
        end
    end
    return c
end

function gene_reaction_rules(model::JSONModel)
    grrs = String[]
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            for rxn in model.m[k]
                push!(grrs, rxn["gene_reaction_rules"]) # assume only key
            end
        end
    end 
end
