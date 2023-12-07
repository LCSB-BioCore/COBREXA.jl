
#TODO: at this point, consider renaming the whole thing to "settings"

"""
$(TYPEDSIGNATURES)

Change the objective sense of optimization. Possible arguments are
`JuMP.MAX_SENSE` and `JuMP.MIN_SENSE`.
"""
set_objective_sense(objective_sense) =
    opt_model -> J.set_objective_sense(opt_model, objective_sense)

export set_objective_sense

"""
$(TYPEDSIGNATURES)

Change the JuMP optimizer used to run the optimization.
"""
set_optimizer(optimizer) = opt_model -> J.set_optimizer(opt_model, optimizer)

export set_optimizer

"""
$(TYPEDSIGNATURES)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
set_optimizer_attribute(attribute_key, value) =
    opt_model -> J.set_optimizer_attribute(opt_model, attribute_key, value)

export set_optimizer_attribute

"""
    silence

Modification that disable all output from the JuMP optimizer (shortcut for
`set_silent` from JuMP).
"""
silence(opt_model) = J.set_silent(opt_model)

export silence
