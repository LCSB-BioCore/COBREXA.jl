"""
ReactionFlux associates a rxn with a flux.
"""
mutable struct ReactionFlux
    rxn :: Reaction
    flux :: Float64
end

"""
ReactionFluxes holds references to the objective, reactions and fluxes from analysis (e.g. FBA)
"""
mutable struct ReactionFluxes
    objective_id :: String 
    objective :: Float64
    rxnfluxes :: Array{ReactionFlux, 1}
end

"""
getindex(reactionfluxes, rxn)

Return the index of rxn in reactionfluxes and -1 if it is not found.
Note, this is slightly different from the normal getindex function.
"""
function Base.getindex(rfs::ReactionFluxes, rxn::Reaction)
    for i in eachindex(rfs.rxnfluxes)
        if rxn.id == rfs.rxnfluxes[i].rxn.id
            return i
        end
    end
    return -1
end

"""
Pretty printing of ReactionFluxes objects.
"""
function Base.show(io::IO, rfs::ReactionFluxes)
    inds = sortperm(abs.([rf.flux for rf in rfs.rxnfluxes]), rev=true) # max abs fluxes 
    println(io, "Optimum for $(rfs.objective_id) = ", round(rfs.objective, digits=4))
    counter = 0
    for i in inds
        if startswith(rfs.rxnfluxes[i].rxn.id, "EX_")
            println(io, rfs.rxnfluxes[i].rxn.name, " = ", round(rfs.rxnfluxes[i].flux, digits=4), " mmol/gDW/h")
            counter += 1
        end
        if counter > 10 # only display top 10
            break
        end
    end
end

"""
cbmodel = CBM(model::Model)

Initialize a constraint based model. Creates a model that satisfies the mass balance
and flux constraints but no objective or optimizer is set. Returns the JuMP model.
"""
function CBM(model::Model)
    coremodel = CoreModel(model) # Construct S, b, lbs, ubs from model
    cbmodel = JuMP.Model()
    nvars = size(coremodel.S, 2) # number of variables in model
    v = @variable(cbmodel, v[1:nvars]) # flux variables
    mb = @constraint(cbmodel, mb, coremodel.S*v .== coremodel.b) # mass balance
    lbs = @constraint(cbmodel, lbs, coremodel.lbs .<= v) # lower bounds
    ubs = @constraint(cbmodel, ubs, v .<= coremodel.ubs) # upper bounds
    return cbmodel, v, mb, ubs, lbs
end

"""
reactionfluxes = fba(model::Model, objective_rxns; weights = Float64[], optimizer="gubori")

Run flux balance analysis on the model using objective_rxn(s) and optionally specifying their weights.
Optimiser can also be set here. Uses the constraints implied by the model object. Objective is set separately.
"""
function fba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")
    cbm, _, _, _, _ = CBM(model) # get the base constraint based model

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
        objective_id = length(objective_indices) == 1 ? model.rxns[objective_indices[1]].id : "multiple reactions"
        arr= Array{ReactionFlux, 1}()
        for (i, rxn) in enumerate(model.rxns)
            push!(arr, ReactionFlux(rxn, value(v[i])))
        end
        return ReactionFluxes(objective_id, objective_value(cbm), arr)        
    else
        return ReactionFluxes("Optimization issues", 0.0, Array{ReactionFlux, 1}())
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
    cbm, _, _, _, _ = CBM(model) # get the base constraint based model

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
        arr = Array{ReactionFlux, 1}()
        for (i, rxn) in enumerate(model.rxns)
            push!(arr, ReactionFlux(rxn, value(v[i])))
        end
        return ReactionFluxes("Σᵢ||vᵢ||", objective_value(cbm), arr)        
    else
        return ReactionFluxes("Optimization issues", 0.0, Array{ReactionFlux, 1}())
    end
end

"""
atom_balance_dict = atom_exchange(reactionfluxes)

Return the composition of atoms consumed or produced by the model according to reactionfluxes.
"""
function atom_exchange(rfs::ReactionFluxes)
    # find exchange reactions
    ex_inds = [i for i in eachindex(rfs.rxnfluxes) if startswith(rfs.rxnfluxes[i].rxn.id, "EX_")]
    
    atom_balance = Dict{String, Float64}()
    for ex_ind in ex_inds
        for (met, w) in rfs.rxnfluxes[ex_ind].rxn.metabolites
            for (atom, stoich) in getatoms(met)
                atom_balance[atom] = get(atom_balance, atom, 0.0) + stoich*w*value(rfs.rxnfluxes[ex_ind].flux)
            end
        end
    end
    return atom_balance
end

function map_fluxes(v, model)
    rxndict = Dict{String, Float64}()
    for i in eachindex(model.rxns)
        rxndict[model.rxns[i].id] = value(v[i])
    end
    return rxndict
end
