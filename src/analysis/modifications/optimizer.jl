
"""
    change_sense(objective_sense)

Change the objective sense of optimization.
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.

If you want to change the objective and sense at the same time, use
[`change_objective`](@ref) instead to do both at once.
"""
function change_sense(objective_sense)
    (_, opt_model) -> set_objective_sense(opt_model, objective_sense)
end

"""
    change_optimizer(optimizer)

Change the JuMP optimizer used to run the optimization.

This may be used to try different approaches for reaching the optimum, and in
problems that may require different optimizers for different parts, such as the
[`parsimonious_flux_balance_analysis`](@ref).
"""
function change_optimizer(optimizer)
    (_, opt_model) -> set_optimizer(opt_model, optimizer)
end

"""
    change_optimizer_attribute(attribute_key, value)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
function change_optimizer_attribute(attribute_key, value)
    (_, opt_model) ->
        set_optimizer_attribute(opt_model, attribute_key, value)
end
