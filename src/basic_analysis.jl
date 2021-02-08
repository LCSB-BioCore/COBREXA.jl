"""
cbmodel = initCBM(coremodel :: CoreModel; optimizer="Gurobi")

Initialize a constraint based model. Creates a model that satisfies the mass balance
and flux constraints but no objective is set. Return references to these objects for 
simple modification if necessary.
"""
function initCBM(coremodel :: CoreModel; optimizer="gurobi")
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
    
    return cbmodel, v, massbalance, fluxlbs, fluxubs
end

"""
fba(model::Model, objfunctionrxn; optimizer="gubori")

Run flux balance analysis on the Model
"""
function fba(model :: Model, objective_rxn; optimizer="gurobi")
    coremodel = CoreModel(model)
    cbmodel, v, massbalance, fluxlbs, fluxubs = initCBM(coremodel, optimizer=optimizer)
    
    objective_index = model[objective_rxn]
    @objective(cbmodel, Max, v[objective_index])
    optimize!(cbmodel)
    cto.verbose && @info "FBA status: $(termination_status(cbmodel.cbmodel))"
    return v
end


# μ = objective_value(model) 

# # Solve pFBA problem
# @constraint(model, 0.999*μ <= v[obj_ind] <= μ) # set biomass function to FBA solution
# @objective(model, Min, sum(dot(v,v)))
# optimize!(model)
# println("pFBA status: ", termination_status(model))

# return v, μ
