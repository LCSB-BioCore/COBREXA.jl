"""
cbmodel = initCBM(model :: Model; optimizer="Gurobi")

Initialize a constraint based model. Creates a model that satisfies the mass balance
and flux constraints but no objective is set. Return references to these objects for 
simple modification if necessary.
"""
function initCBM(model :: Model; optimizer="gurobi")
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
        verbose && @warn "Optimizer not yet directly supported, however, you can set it yourself with `set_optimizer(cbmodel, OPTIMIZER)`.\nSee JuMP's documentation."
    end
    
    nvars = size(model.coremodel.S, 2) # number of variables in model
    
    v = @variable(cbmodel, v[1:nvars]) # flux variables
    @constraint(cbmodel, massbalance, model.coremodel.S*v .== model.coremodel.b) # mass balance constraints
    @constraint(cbmodel, fluxlbs, model.coremodel.lbs .<= v)
    @constraint(cbmodel, fluxubs, v .<= model.coremodel.ubs)

    cbmodelout = (cbmodel=cbmodel, v=v, massbalance=massbalance, fluxlbs=fluxlbs, fluxubs=fluxubs) # named tuple for convenience
    
    return cbmodelout
end

function fba(cbmodel, objective_index)
    @objective(cbmodel.cbmodel, Max, cbmodel.v[objective_index])
    optimize!(cbmodel.cbmodel)
    verbose && @info "FBA status: $(termination_status(cbmodel.cbmodel))"
end

"""
index = getindex(model, rxn)

Get the index of rxn in model
"""
function getindex(model::Model, rxn::Reaction)
    for i in eachindex(model.rxns)
        if model.rxns[i].id == rxn.id
            return i
        end
    end
    return -1
end


"""
"""
function findrxn(model::Model, rxnid::String)
    for i in eachindex(model.rxns)
        if model.rxns[i].id == rxnid
            return model.rxns[i]
        end
    end
    return nothing
end

# μ = objective_value(model) 

# # Solve pFBA problem
# @constraint(model, 0.999*μ <= v[obj_ind] <= μ) # set biomass function to FBA solution
# @objective(model, Min, sum(dot(v,v)))
# optimize!(model)
# println("pFBA status: ", termination_status(model))

# return v, μ
