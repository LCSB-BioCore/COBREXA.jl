function fba(model :: Model)

end

function fba(model :: CoreModel, objective_index)

end

nvars = size(rawmodel["S"], 2) # number of variables in model
model = Model(Gurobi.Optimizer) # model
set_optimizer_attribute(model, "OutputFlag", 0) # quiet

v = @variable(model, v[1:nvars]) # flux variables
@constraint(model, massbalance, rawmodel["S"]*v .== rawmodel["b"]) # mass balance constraints

constraint_subset_inds = getrxninds(rawmodel, keys(constraint_subset))
default_constraints = filter(x -> !(x in values(constraint_subset_inds)), 1:length(rawmodel["rxns"]))
@constraint(model, fluxbounds, rawmodel["lb"][default_constraints] .<= v[default_constraints] .<= rawmodel["ub"][default_constraints]) # flux bounds


for (k, cv) in constraint_subset
    @constraint(model, cv[1] <= v[constraint_subset_inds[k]] <= cv[2])
end

return model, v

model, v = mkCBM(rawmodel, constraint_subset)
    
# Solve FBA problem
obj_ind = getrxnind(rawmodel, objective_id)
@objective(model, Max, v[obj_ind]) # objective biomass maximization
optimize!(model)
println("FBA status: ", termination_status(model))

μ = objective_value(model) 

# Solve pFBA problem
@constraint(model, 0.999*μ <= v[obj_ind] <= μ) # set biomass function to FBA solution
@objective(model, Min, sum(dot(v,v)))
optimize!(model)
println("pFBA status: ", termination_status(model))

return v, μ
