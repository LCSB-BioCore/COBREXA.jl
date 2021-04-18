struct JSONModel <: MetabolicModel
    m::Dict{String,Any}
end

# Generic interface
# Unfortunately this model type does not have standardized field names, hence the need to look for valid fieldnames.
# The keys used to look for valid fieldnames is in `constants`.
# However, assume that `id` is used to index reactions, metabolites and genes.
# Also assume that `metabolites` is used to access the dict containing the reaction equation for a reaction.

function reactions(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in possible_rxn_keys
        if haskey(model.m, k)
            return [string(r["id"]) for r in model.m[k][:]]
        end
    end
    @warn "No reactions found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_reactions(model::JSONModel)::Int = length(reactions(model))

function metabolites(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in possible_met_keys
        if haskey(model.m, k)
            return [string(m["id"]) for m in model.m[k][:]]
        end
    end
    @warn "No metabolites found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_metabolites(model::JSONModel)::Int = length(metabolites(model))

function genes(model::JSONModel)::Union{Nothing, Vector{String}}
    for k in possible_gene_keys
        if haskey(model.m, k)
            return [string(g["id"]) for g in model.m[k][:]]
        end
    end
    @warn "No genes found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_genes(model::JSONModel)::Int = length(genes(model))

function stoichiometry(model::JSONModel)::Union{Nothing, SparseMat}
    rxn_ids = reactions(model)
    rxn_key = _constants.possible_rxn_keys[[haskey(model.m, x) for x in _constants.possible_rxn_keys]][1] # get the rxn key used
    met_ids = metabolites(model)
    S = SparseArrays.spzeros(length(mets), length(rxn))
    for (i, rxn_id) in rxn_ids
        rxn_dict = model.m[rxn_key][i]["metabolites"]
        for (met_id, coeff) in rxn_dict
            j = findfirst(x -> x == met_id, met_ids) # row
            isnothing(j) ?
            (@error "S matrix construction error: $(met_id) not defined."; continue) : (return nothing)
            S[j, i] = coeff
        end
    end
    return S
end

function lower_bounds(model::JSONModel)
    ks = ("lbs", "lb", "lowerbounds", "lower_bounds")
    for k in ks
        if haskey(model.m, k)
            return [float(r) for r in model.m[k][:]]
        end
    end
    @warn "No lower bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function upper_bounds(model::JSONModel)
    ks = ("ubs", "ub", "upperbounds", "upper_bounds")
    for k in ks
        if haskey(model.m, k)
            return [float(r) for r in model.m[k][:]]
        end
    end
    @warn "No upper bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function bounds(model::JSONModel)
    return lower_bounds(model), upper_bounds(model)
end

function balance(model::JSONModel)
    ks = ("b")
    for k in ks
        if haskey(model.m, k)
            return sparse([float(x) for x in model.m[k][:]])
        end
    end
    return spzeros(length(model.metabolites))
end

function objective(model::JSONModel)
    ks = ("c")
    for k in ks
        if haskey(model.m, k)
            return sparse([float(x) for x in model.m[k][:]])
        end
    end
    @warn "No objective vector found. Perhaps the an exotic field name is used by the model?"
    return nothing
end
