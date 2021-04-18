"""
    struct MATModel

Model struct used when importing models saved in `.mat` format.
All fields in model are imported (i.e. no data loss occurs).
However, not all the fields are used by analysis functions.
"""
struct MATModel <: MetabolicModel
    id::String
    m::Dict{String,Any}
end

# Generic interface
# Unfortunately this model type does not have standardized field names, hence the need to look for valid fieldnames.
# The keys used to look for valid fieldnames is in `constants`.

function reactions(model::MATModel)::Union{Nothing,Vector{String}}
    for k in _constants.possible_rxn_keys
        if haskey(model.m, k)
            return [string(r) for r in model.m[k][:]] # sometimes stored as a matrix, this ensure that it is a string vector
        end
    end
    @warn "No reactions found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_reactions(model::MATModel)::Int = length(reactions(model))

function metabolites(model::MATModel)::Union{Nothing,Vector{String}}
    for k in _constants.possible_met_keys
        if haskey(model.m, k)
            return [string(r) for r in model.m[k][:]] # sometimes stored as a matrix, this ensure that it is a string vector
        end
    end
    @warn "No metabolites found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_metabolites(model::MATModel)::Int = length(metabolites(model))

function genes(model::MATModel)::Union{Nothing,Vector{String}}
    for k in _constants.possible_gene_keys
        if haskey(model.m, k)
            return [string(r) for r in model.m[k][:]] # sometimes stored as a matrix, this ensure that it is a string vector
        end
    end
    @warn "No genes found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

n_genes(model::MATModel)::Int = length(genes(model))

function stoichiometry(model::MATModel)::Union{Nothing,SparseMat}
    for k in _constants.possible_stoich_matrix_keys
        if haskey(model.m, k)
            return sparse(model.m[k])
        end
    end
    @warn "No stoichiometric matrix found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function lower_bounds(model::MATModel)
    for k in _constants.possible_lower_bound_keys
        if haskey(model.m, k)
            return [float(r) for r in model.m[k][:]]
        end
    end
    @warn "No lower bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function upper_bounds(model::MATModel)
    for k in _constants.possible_upper_bound_keys
        if haskey(model.m, k)
            return [float(r) for r in model.m[k][:]]
        end
    end
    @warn "No upper bounds found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function bounds(model::MATModel)
    return lower_bounds(model), upper_bounds(model)
end

function balance(model::MATModel)
    for k in _constants.possible_balance_keys
        if haskey(model.m, k)
            return sparse([float(x) for x in model.m[k][:]])
        end
    end
    return spzeros(n_metabolites(model))
end

function objective(model::MATModel)
    for k in _constants.possible_objective_keys
        if haskey(model.m, k)
            return sparse([float(x) for x in model.m[k][:]])
        end
    end
    @warn "No objective vector found. Perhaps the an exotic field name is used by the model?"
    return spzeros(n_reactions(model))
end

function id(model::MATModel)
    return model.id
end

function gene_reaction_rules(model::MATModel)
    for k in _constants.possible_grr_keys
        if haskey(model.m, k)
            return [string(x) for x in model.m[k][:]]
        end
    end
    @warn "No gene reaction rules found. Perhaps the an exotic field name is used by the model?"
    return nothing
end

function build_reactions(model::MATModel)

end

function build_genes(model::MATModel)

end

function build_metabolites(model::MATModel)

end
