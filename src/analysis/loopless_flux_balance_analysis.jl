function _get_boundary_reactions(model::StandardModel)
    return [i for i in model.reactions if length(i.metabolites) == 1]
end

function _set_solver_attibutes!(opt_model, solver_attributes)
    if !isempty(solver_attributes) # set other attributes
        for (k, val) in solver_attributes
            set_optimizer_attribute(opt_model, k, val)
        end
    end
end

function _set_additional_constraints!(model, opt_model, constraints)
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, opt_model; lb=con[1], ub=con[2])
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

function _constrain_objective_value!(opt_model, opt_weights, objective_indices)
    λ = objective_value(opt_model)
    v = opt_model[:x]
    @constraint(
        opt_model,
        λ <= sum(opt_weights[i] * v[i] for i in objective_indices) <= λ
    )
end

function _add_cycle_free!(model, opt_model, objective_indices, fluxes)
    v = opt_model[:x]
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
        set_bound(i, opt_model; ub=ub, lb=lb)
    end    
    @objective(opt_model, Min, dot(min_objectives, v))
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
    model::M,
    opt_model,
    opt_weights,
    objective_indices,
    fluxes,
)::Union{Dict{String,Float64},Nothing} where {M <: MetabolicModel}
    _constrain_objective_value!(opt_model, opt_weights, objective_indices)
    _add_cycle_free!(model, opt_model, objective_indices, fluxes)
    optimize!(opt_model)
    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing
    return Dict(zip(reactions(model), value.(opt_model[:x])))
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
function loopless_flux_balance_analysis(
    model::M,
    optimizer;
    objective_func::Union{Reaction,Array{Reaction,1}},
    weights=[],
    modifications=[(model, opt_model) -> nothing],

)::Union{Dict{String,Float64},Nothing} where {M <: MetabolicModel}
    # FBA PART
    opt_model = make_optimization_model(model, optimizer)

    # Is there a way to get the objective indices if we set the 
    # objective via `change_objective`? then we could recycle the flux_balance_analysis
    # method in a better way
    objective_indices = _get_reaction_indices(model, objective_func)
    opt_weights = _get_objective_weights(model, objective_indices, weights)
    change_objective(objective_func, weights=weights)(model, opt_model)

    # support for multiple modification, fallback to single one
    if typeof(modifications) <: AbstractVector
        for mod in modifications
            mod(model, opt_model)
        end
    else
        modifications(model, opt_model)
    end

    optimize!(opt_model)
    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing
    
    fluxes = Dict(zip(reactions(model), value.(opt_model[:x])))
    return loopless_solution(
        model,
        opt_model,
        opt_weights,
        objective_indices,
        fluxes,
    )
end



