"""
cbmodel, v, mb, ubs, lbs = CBM(model::Model)

Initialize a constraint based model. Creates a model that satisfies the mass balance
and flux constraints but no objective or optimizer is set. Returns the JuMP model.
This is useful if you want to write your own optimization problem.

cbmodel is the JuMP model. v are the fluxes, mb is S*v == 0, and lbs <= v <= ubs.
"""
function build_cbm(model::Model)
    S, b, ubs, lbs = get_core_model(model) # Construct S, b, lbs, ubs from model
    cbmodel = JuMP.Model()
    nvars = size(S, 2) # number of variables in model
    v = @variable(cbmodel, v[1:nvars]) # flux variables
    mb = @constraint(cbmodel, mb, S*v .== b) # mass balance
    lbs = @constraint(cbmodel, lbs, lbs .<= v) # lower bounds
    ubs = @constraint(cbmodel, ubs, v .<= ubs) # upper bounds
    return cbmodel, v, mb, ubs, lbs
end

"""
reactionfluxes = fba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights = Float64[], optimizer="gubori")

Run flux balance analysis (FBA) on the model using objective_rxn(s) and optionally specifying their weights (empty weights mean equal weighting per reaction).
Optimiser can also be set here. Uses the constraints implied by the model object.
"""
function fba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")
    cbm, _, _, _, _ = build_cbm(model) # get the base constraint based model

    # ensure that an array of objective indices are fed in
    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end
    
    # ensure corrects weights are given to objective
    if isempty(weights)
        weights = ones(length(objective_indices))
    end
    
    # update the objective function tracker
    wcounter = 1
    for i in eachindex(model.rxns)
        if i in objective_indices
            model.rxns[i].objective_coefficient = weights[wcounter]
            wcounter += 1
        else
            model.rxns[i].objective_coefficient = 0.0
        end
    end

    # set optimizer
    if optimizer == "glpk"
        set_optimizer(cbm, GLPK.Optimizer)
        set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF) # quiet
    elseif optimizer == "gurobi"
        set_optimizer(cbm, Gurobi.Optimizer)
        set_optimizer_attribute(cbm, "OutputFlag", 0) # quiet
    elseif optimizer == "tulip"
        set_optimizer(cbm, Tulip.Optimizer) # quiet by default
    elseif optimizer == "ipopt"
        set_optimizer(cbm, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0) # quiet    
    else
        @warn "Optimizer not yet directly supported, however, you can set it yourself with `set_optimizer(cbmodel, OPTIMIZER)`.\nSee JuMP's documentation."
    end
    
    v = all_variables(cbm)
    @objective(cbm, Max, sum(v[i] for i in objective_indices))
    optimize!(cbm)    

    status = termination_status(cbm) == MOI.OPTIMAL
    
    if status
        return map_fluxes(v, model.rxns)      
    else
        return Dict{String, Float64}()
    end
end

"""
reactionfluxes = pfba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")

Run parsimonious flux balance analysis on the model using objective_rxn(s) as the objective for the initial FBA problem. 
Optionally, this/these objectives can be weighted by weights.
Optimiser can also be set here. Uses the constraints implied by the model object.
"""
function pfba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")
    ## FBA ################################################
    cbm, _, _, _, _ = build_cbm(model) # get the base constraint based model

    # ensure that an array of objective indices are fed in
    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end
    
    # ensure corrects weights are given to objective
    if isempty(weights)
        weights = ones(length(objective_indices))
    end
    
    # update the objective function tracker
    wcounter = 1
    for i in eachindex(model.rxns)
        if i in objective_indices
            model.rxns[i].objective_coefficient = weights[wcounter]
            wcounter += 1
        else
            model.rxns[i].objective_coefficient = 0.0
        end
    end

    # set optimizer
    if optimizer == "glpk"
        set_optimizer(cbm, GLPK.Optimizer)
        set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF) # quiet
    elseif optimizer == "gurobi"
        set_optimizer(cbm, Gurobi.Optimizer)
        set_optimizer_attribute(cbm, "OutputFlag", 0) # quiet
    elseif optimizer == "tulip"
        set_optimizer(cbm, Tulip.Optimizer) # quiet by default
    elseif optimizer == "ipopt"
        set_optimizer(cbm, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0) # quiet    
    else
        @warn "Optimizer not yet directly supported, however, you can set it yourself with `set_optimizer(cbmodel, OPTIMIZER)`.\nSee JuMP's documentation."
    end
    
    v = all_variables(cbm)
    @objective(cbm, Max, sum(v[i] for i in objective_indices))
    optimize!(cbm)    

    fba_status = termination_status(cbm) == MOI.OPTIMAL

    ## pFBA ###############################################
    λ = objective_value(cbm)

    @constraint(cbm, pfbacon, 0.999999*λ <= sum(v[i] for i in objective_indices) <= λ) # constrain model - 0.9999 should be close enough?
    @objective(cbm, Min, sum(dot(v, v)))
    optimize!(cbm)

    if termination_status(cbm) != MOI.OPTIMAL # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.99999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end
    if termination_status(cbm) != MOI.OPTIMAL # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.9999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end
    if termination_status(cbm) != MOI.OPTIMAL # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end

    pfba_status = termination_status(cbm) == MOI.OPTIMAL
    
    if fba_status && pfba_status
        return map_fluxes(v, model.rxns)      
    else
        @warn "Optimization issues occurred."
        return Dict{String, Float64}() # something went wrong
    end
end

"""
rxn_flux_dict = map_fluxes(v, model::Model)

Map fluxes from an optimization problem (v) to rxns in a model. v can be a JuMP object (fluxes) or an array of Float64 fluxes.
Assumes they are in order, which they should be since they are constructed from model.
"""
function map_fluxes(v::Array{VariableRef,1}, model::Model)
    rxndict = Dict{String, Float64}()
    for i in eachindex(model.rxns)
        rxndict[model.rxns[i].id] = value(v[i])
    end
    return rxndict
end

function map_fluxes(v::Array{Float64,1}, model::Model)
    rxndict = Dict{String, Float64}()
    for i in eachindex(model.rxns)
        rxndict[model.rxns[i].id] = v[i]
    end
    return rxndict
end

"""
setbound(index, ubconstaintref, lbconstaintref; ub=1000, lb=-1000)

Helper function to set the bounds of variables.
The JuMP set_normalized_rhs function is a little confusing, so this function simplifies setting
constraints. Just supply the constraint index and the bound objects and they will be set appropriately.
"""
function set_bound(vind, ubs, lbs; ub=1000, lb=-1000)
    if lb <= 0 
        set_normalized_rhs(lbs[vind], abs(lb))
    else
        set_normalized_rhs(lbs[vind], -abs(lb))
    end
    set_normalized_rhs(ubs[vind], ub)
end

"""
get_exchanges(rxndict::Dict{String, Float64}; topN=8, ignorebound=1000)

Display the topN producing and consuming exchange fluxes. Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
"""
function display_exchange_reactions(rxndict::Dict{String, Float64}; topN=8, ignorebound=1000)
    fluxes = Float64[]
    rxns = String[]
    for (k, v) in rxndict
        if startswith(k, "EX_") && abs(v) < ignorebound
            push!(rxns, k)
            push!(fluxes, v)
        end
    end
    inds_prod = sortperm(fluxes, rev=true)
    inds_cons = sortperm(fluxes)

    println("Consuming fluxes:")
    for i in 1:topN
        if rxndict[rxns[inds_cons[i]]] > 0
            continue
        end
        println(rxns[inds_cons[i]], " = ", round(rxndict[rxns[inds_cons[i]]], digits=4))
    end
    println("Producing fluxes:")
    for i in 1:topN
        if rxndict[rxns[inds_cons[i]]] < 0
            continue
        end
        println(rxns[inds_prod[i]], " = ", round(rxndict[rxns[inds_prod[i]]], digits=4))
    end
end

"""
"""
function display_metabolite_fluxes(fluxdict::Dict{String, Float64}, model::Model)

end
