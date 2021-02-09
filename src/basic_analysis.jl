"""
ReactionFluxes holds references to reactions, fluxes from analysis (e.g. FBA)
"""
mutable struct ReactionFluxes
    rxns :: Array{Reaction, 1} # each metabolite has an associated concentration
    objective_id :: String 
    objective :: Float64
    fluxes :: Array{Float64, 1} # order will match because CoreModel is constructed from rxns in order
end

"""
updatefluxes!(reactionData, newfluxes)

Sets the fluxes field of reactionData to the new fluxes
"""
function update!(rd::ReactionFluxes, objective_id::String, objective::Float64,  fluxes::Array{Float64, 1})
    rd.objective_id = objective_id
    rd.objective = objective
    for i in eachindex(rd)
        rd.fluxes[i] = fluxes[i]
    end
end

function Base.show(io::IO, rd::ReactionFluxes)
    afluxes = abs.(rd.fluxes)
    if maximum(afluxes) ≈ 0.0
        println(io, "No solution found.")
    else
        inds = sortperm(afluxes, rev=true) # display maximum absolute fluxes 
        println(io, "Optimum for $(rd.objective_id) = ", round(rd.objective, digits=4))
        counter = 0
        for i in inds
            if startswith(rd.rxns[i].id, "EX_")
                println(io, rd.rxns[i].name, " = ", round(rd.fluxes[i], digits=4), " mmol/gDW/h")
                counter += 1
            end
            if counter > 10
                break
            end
        end
    end
end

"""
cbmodel = initCBM(coremodel :: CoreModel; optimizer="Gurobi")

Initialize a constraint based model. Creates a model that satisfies the mass balance
and flux constraints but no objective is set. Return references to these objects for 
simple modification if necessary.
"""
function initCBM(coremodel::CoreModel; optimizer="gurobi")
    cbmodel = JuMP.Model()

    # set optimizer
    if optimizer == "glpk"
        set_optimizer(cbmodel, GLPK.Optimizer)
        set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF) # quiet
    elseif optimizer == "gurobi"
        set_optimizer(cbmodel, Gurobi.Optimizer)
        set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet
    elseif optimizer == "tulip"
        set_optimizer(cbmodel, Tulip.Optimizer) # quiet by default
    elseif optimizer == "ipopt"
        set_optimizer(cbmodel, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0) # quiet    
    else
        cto.verbose && @warn "Optimizer not yet directly supported, however, you can set it yourself with `set_optimizer(cbmodel, OPTIMIZER)`.\nSee JuMP's documentation."
    end
    
    nvars = size(coremodel.S, 2) # number of variables in model
    
    v = @variable(cbmodel, v[1:nvars]) # flux variables
    @constraint(cbmodel, massbalance, coremodel.S*v .== coremodel.b) # mass balance constraints
    @constraint(cbmodel, fluxlbs, coremodel.lbs .<= v)
    @constraint(cbmodel, fluxubs, v .<= coremodel.ubs)
    
    return cbmodel
end

"""
solvefba(model::Model, objective_rxns; weights = Float64[], optimizer="gubori")

Run flux balance analysis on the model using objective_rxn(s) and optionally specifying their weights.
Optimiser can also be set here.
"""
function solvefirstfba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")

    coremodel = CoreModel(model) # Construct S, b, lbs, ubs from model
    cbm = initCBM(coremodel, optimizer=optimizer) # create the optimization problem

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
    v = all_variables(cbm)
    @objective(cbm, Max, sum(v[i] for i in objective_indices))
    optimize!(cbm)    

    return cbm
end

function fba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")
    
    cbm = solvefirstfba(model, objective_rxns; weights=weights, optimizer=optimizer)
    
    v = all_variables(cbm)

    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end

    status = termination_status(cbm) == MOI.OPTIMAL
    
    if status
        objective_id = length(objective_indices) == 1 ? model.rxns[objective_indices[1]].id : "multiple reactions"    
        return ReactionFluxes(model.rxns, objective_id, objective_value(cbm), [value(v[i]) for i in eachindex(v)])        
    else
        return ReactionFluxes(model.rxns, "", 0.0, zeros(length(v)))
    end
end

"""
solobj = pfba()

Not efficient.
"""
function pfba(model::Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}; weights=Float64[], optimizer="gurobi")
    
    cbm = solvefirstfba(model, objective_rxns; weights=weights, optimizer=optimizer)
    v = all_variables(cbm)

    statusfba = termination_status(cbm) == MOI.OPTIMAL
    
    # ensure that an array of objective indices are fed in
    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end
    
    λ = objective_value(cbm)

    @constraint(cbm, pfbacon, 0.9999*λ <= sum(v[i] for i in objective_indices) <= λ) # constrain model - 0.9999 should be close enough?
    @objective(cbm, Min, sum(dot(v, v)))
    optimize!(cbm)

    if termination_status(cbm) != MOI.OPTIMAL # try to relax bound
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, pfbacon, 0.999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end

    statuspfba = termination_status(cbm) == MOI.OPTIMAL
    
    if statusfba && statuspfba
        return ReactionFluxes(model.rxns, "Σ||v - ̄v||", objective_value(cbm), [value(v[i]) for i in eachindex(v)])        
    else
        return ReactionFluxes(model.rxns, "", 0.0, zeros(length(v)))
    end
end

"""
Perform pFBA and calculate the atom balance across the model
"""
function atom_exchange(rfs::ReactionFluxes)
    # find exchange reactions
    ex_inds = [i for i in eachindex(rfs.rxns) if startswith(rfs.rxns[i].id, "EX_")]
    
    atom_balance = Dict{String, Float64}()
    for ex_ind in ex_inds
        for (met, w) in rfs.rxns[ex_ind].metabolites
            for (atom, stoich) in getatoms(met)
                atom_balance[atom] = get(atom_balance, atom, 0.0) + stoich*w*value(rfs.fluxes[ex_ind])
            end
        end
    end
    return atom_balance
end
