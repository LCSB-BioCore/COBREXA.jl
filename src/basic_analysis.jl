struct Solution
    isoptimal :: Bool
    objid :: String
    obj :: Float64
    fluxids :: Array{String, 1}
    fluxnames :: Array{String, 1}
    fluxes :: Array{Float64, 1}
end

function Base.show(io::IO, s::Solution)
    afluxes = abs.(s.fluxes)
    inds = sortperm(afluxes, rev=true)
    if s.isoptimal
        println("Optimum for $(s.objid) = ", round(s.obj, digits=4))
        counter = 0
        for i in inds
            if startswith(s.fluxids[i], "EX_")
                println(s.fluxnames[i], " = ", round(s.fluxes[i], digits=4), " mmol/gDW/h")
                counter += 1
            end
            if counter > 10
                break
            end
        end
    else
        println("Optimization issues, try again.")
    end
end

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

Run flux balance analysis on the Model. Inefficient implementation.
"""
function fba(model :: Model, objective_rxn; optimizer="gurobi")
    objective_index = model[objective_rxn]
    cbmodel, v = dofba(model, objective_rxn; optimizer)

    status = termination_status(cbmodel) == MOI.OPTIMAL
    
    solobj = Solution(status, objective_rxn.id, value(v[objective_index]), [rxn.id for rxn in model.rxns], [rxn.name for rxn in model.rxns], [value(v[i]) for i in eachindex(v)])

    return solobj
end

function dofba(model :: Model, objective_rxn; optimizer="gurobi")
    coremodel = CoreModel(model)
    cbmodel, v, massbalance, fluxlbs, fluxubs = initCBM(coremodel, optimizer=optimizer)
    
    objective_index = model[objective_rxn]
    @objective(cbmodel, Max, v[objective_index])
    optimize!(cbmodel)
    cto.verbose && @info "FBA status: $(termination_status(cbmodel))"

    return cbmodel, v
end


"""
solobj = pfba()

Not efficient.
"""
function pfba(model::Model, objective_rxn; optimizer="gurobi")
    objective_index = model[objective_rxn]
    cbmodel, v = dopfba(model, objective_rxn; optimizer)

    status = termination_status(cbmodel) == MOI.OPTIMAL
    
    solobj = Solution(status, objective_rxn.id, value(v[objective_index]), [rxn.id for rxn in model.rxns], [rxn.name for rxn in model.rxns], [value(v[i]) for i in eachindex(v)])

    return solobj
end

function dopfba(model::Model, objective_rxn; optimizer="gurobi")

    objective_index = model[objective_rxn]
    cbmodel, v = dofba(model, objective_rxn; optimizer)
    λ = value(v[objective_index])
    @constraint(cbmodel, 0.9999*λ <= v[objective_index] <= λ) # constrain model - 0.9999 should be close enough?
    @objective(cbmodel, Min, sum(dot(v, v)))
    optimize!(cbmodel)
    cto.verbose && @info "pFBA status: $(termination_status(cbmodel))"
    return cbmodel, v
end



