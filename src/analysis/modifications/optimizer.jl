
"""
$(TYPEDSIGNATURES)

Change the objective sense of optimization. Possible arguments are
`JuMP.MAX_SENSE` and `JuMP.MIN_SENSE`.
"""
modify_sense(objective_sense) =
    (_, opt_model) -> set_objective_sense(opt_model, objective_sense)

"""
$(TYPEDSIGNATURES)

Change the JuMP optimizer used to run the optimization.
"""
modify_optimizer(optimizer) = (_, opt_model) -> J.set_optimizer(opt_model, optimizer)

"""
$(TYPEDSIGNATURES)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
modify_optimizer_attribute(attribute_key, value) =
    (_, opt_model) -> J.set_optimizer_attribute(opt_model, attribute_key, value)

"""
    silence

Modification that disable all output from the JuMP optimizer (shortcut for
`set_silent` from JuMP).
"""
const silence = (_, opt_model) -> J.set_silent(opt_model)

export modify_sense, modify_optimizer, modify_optimizer_attribute, silence
