
"""
    change_sense(objective_sense)

Change the objective sense. 
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.
Note, [`change_objective`](@ref) sets the sense of the objective,
so it doesn't make sense to use this function AND [`change_objective`](@ref) simultaneously. 
"""
function change_sense(objective_sense)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_objective_sense(opt_model, objective_sense)
    end
end

"""
    change_solver(optimizer)

Change the solver (`optimizer`) used to solve the model.
Typically the solver is specified as a required argument in a function. 
However, this function is useful if the problem has multiple subparts that 
require different solvers.

See also: [`parsimonious_flux_balance_analysis`](@ref)
"""
function change_solver(optimizer)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer(opt_model, optimizer)
    end
end

"""
    change_solver_attribute(option_key, option_val)

Change a solver attribute. These attributes are solver specific,
refer the either JuMP or the solver you are using's documentation.
"""
function change_solver_attribute(option_key, option_val)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer_attribute(opt_model, option_key, option_val)
    end
end
