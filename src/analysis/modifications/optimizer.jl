
"""
$(TYPEDSIGNATURES)

Change the objective sense of optimization.
Possible arguments are `MAX_SENSE` and `MIN_SENSE`.

If you want to change the objective and sense at the same time, use
[`change_objective`](@ref) instead to do both at once.
"""
change_sense(objective_sense) =
    (_, opt_model) -> set_objective_sense(opt_model, objective_sense)

"""
$(TYPEDSIGNATURES)

Change the JuMP optimizer used to run the optimization.

This may be used to try different approaches for reaching the optimum, and in
problems that may require different optimizers for different parts, such as the
[`parsimonious_flux_balance_analysis`](@ref).
"""
change_optimizer(optimizer) = (_, opt_model) -> set_optimizer(opt_model, optimizer)

"""
$(TYPEDSIGNATURES)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
change_optimizer_attribute(attribute_key, value) =
    (_, opt_model) -> set_optimizer_attribute(opt_model, attribute_key, value)

"""
    silence

Modification that disable all output from the JuMP optimizer (shortcut for
`set_silent` from JuMP).
"""
const silence = (_, opt_model) -> set_silent(opt_model)
