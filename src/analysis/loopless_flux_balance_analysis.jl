function _get_boundary_reactions(model::CobraModel)
    return [i for i in model.reactions if length(i.metabolites) == 1]
end

function _set_solver_attibutes!(cbm, solver_attributes)
    if !isempty(solver_attributes) # set other attributes
        for (k, val) in solver_attributes
            set_optimizer_attribute(cbm, k, val)
        end
    end
end

function _set_additional_constraints!(model, cbm, constraints)
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, cbm; lb=con[1], ub=con[2])
    end
end

function _get_reaction_indices(model, objective_func)
    if typeof(objective_func) == Reaction
        objective_indices = [model[objective_func]]
    else 
        objective_indices = [model[rxn] for rxn in objective_func]
    end
    return objective_indices
end

function _get_objective_weights(model, objective_indices, weights)
    if isempty(weights)
        weights = ones(length(objective_indices))
    end
    opt_weights = zeros(length(model.reactions))

    # update the objective function tracker
    wcounter = 1
    for i in eachindex(model.reactions)
        if i in objective_indices
            opt_weights[i] = weights[wcounter]
            wcounter += 1
        end
    end
    return opt_weights
end

function _constrain_objective_value!(cbm, opt_weights, v, objective_indices)
    λ = objective_value(cbm)
    @constraint(
        cbm,
        λ <= sum(opt_weights[i] * v[i] for i in objective_indices) <= λ
    )
end

function _add_cycle_free!(model, cbm, objective_indices, v, fluxes)
    boundary_ids = [i.id for i in _get_boundary_reactions(model)]
    min_objectives = zeros(Int64, 1, length(v))
    for (i, reaction) in enumerate(model.reactions)
        id = reaction.id
        if i in objective_indices
            continue
        end
        flux = fluxes[id]
        if id in boundary_ids
            lb = flux
            ub = flux
        elseif flux >= 0
            lb = max(0, reaction.lb)
            ub = max(flux, reaction.ub)
            min_objectives[i] =  1
        else
            lb = min(flux, reaction.lb)
            ub = min(0, reaction.ub)            
            min_objectives[i] = -1
        end
        set_bound(i, cbm; ub=ub, lb=lb)
    end    
    @objective(cbm, Min, dot(min_objectives, v))
end

function _check_status(cbm)
    status = (
        termination_status(cbm) == MOI.OPTIMAL ||
        termination_status(cbm) == MOI.LOCALLY_SOLVED
    )
    if status
        return true
    else
        @warn "Optimization issues occurred."
        return false
    end
end


"""
    loopless_solution(args...)::Union{Dict{String,Float64},Nothing}

Convert an existing FBA solution to a loopless one. Removes as many loops as possible usingthe method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function loopless_solution(
    model::CobraModel,
    cbm,
    opt_weights,
    objective_indices,
    v,
    fluxes,
)
    _constrain_objective_value!(cbm, opt_weights, v, objective_indices)
    _add_cycle_free!(model, cbm, objective_indices, v, fluxes)
    optimize!(cbm)
    if !_check_status(cbm)
        return nothing
    end
    return map_fluxes(v, model)
end

"""
    loopless(args...)::Union{Dict{String,Float64},Nothing}

Run loopless FBA on the model

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function loopless(
    model::CobraModel,
    optimizer,
    objective_func::Union{Reaction,Array{Reaction,1}};
    weights=Float64[],
    solver_attributes=Dict{Any,Any}(),
    constraints=Dict{String,Tuple{Float64,Float64}}(),
    sense=MOI.MAX_SENSE
)
    # FBA PART
    cbm = make_optimization_model(model, optimizer, sense=sense)
    v = cbm[:x]
    _set_solver_attibutes!(cbm, solver_attributes)
    _set_additional_constraints!(model, cbm, constraints)
    objective_indices = _get_reaction_indices(model, objective_func)
    opt_weights = _get_objective_weights(model, objective_indices, weights)
    @objective(cbm, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    optimize!(cbm)
    if !_check_status(cbm)
        return nothing
    end
    fluxes = map_fluxes(v, model)
    return loopless_solution(
        model,
        cbm,
        opt_weights,
        objective_indices,
        v,
        fluxes,
    )
end
