
"""
$(TYPEDSIGNATURES)

Construct a JuMP `Model` that describes the precise constraint system into the
JuMP `Model` created for solving in `optimizer`, with a given optional
`objective` and optimization `sense`.
"""
function optimization_model(
    cs::C.ConstraintTreeElem;
    objective::Union{Nothing,C.Value} = nothing,
    optimizer,
    sense = Maximal,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(cs)])
    isnothing(objective) || J.@objective(model, sense, C.substitute(objective, x))

    # constraints
    function add_constraint(c::C.Constraint)
        if c.bound isa C.EqualTo
            J.@constraint(model, C.substitute(c.value, x) == c.bound.equal_to)
        elseif c.bound isa C.Between
            val = C.substitute(c.value, x)
            isinf(c.bound.lower) || J.@constraint(model, val >= c.bound.lower)
            isinf(c.bound.upper) || J.@constraint(model, val <= c.bound.upper)
        elseif c.bound isa Binary
            anon_bool = J.@variable(model, binary = true)
            J.@constraint(model, C.substitute(c.value, x) == anon_bool)
        end
    end
    function add_constraint(c::C.ConstraintTree)
        add_constraint.(values(c))
    end
    add_constraint(cs)

    return model
end

export optimization_model

"""
$(TYPEDSIGNATURES)

`true` if `opt_model` solved successfully (solution is optimal or
locally optimal). `false` if any other termination status is reached.
"""
is_solved(opt_model::J.Model) =
    J.termination_status(opt_model) in [J.MOI.OPTIMAL, J.MOI.LOCALLY_SOLVED]

export is_solved

"""
    Minimal

Objective sense for finding the minimal value of the objective.

Same as `JuMP.MIN_SENSE`.
"""
const Minimal = J.MIN_SENSE

"""
    Maximal

Objective sense for finding the maximal value of the objective.

Same as `JuMP.MAX_SENSE`.
"""
const Maximal = J.MAX_SENSE

"""
    Maximal

Objective sense for finding the any feasible value of the objective.

Same as `JuMP.FEASIBILITY_SENSE`.
"""
const Feasible = J.FEASIBILITY_SENSE

"""
$(TYPEDSIGNATURES)

Make an JuMP model out of `constraints` using [`optimization_model`](@ref)
(most arguments are forwarded there), then apply the `modifications`, optimize
the model, and return either `nothing` if the optimization failed, or `output`
substituted with the solved values (`output` defaults to `constraints`.

For a "nice" version for simpler finding of metabolic model optima, use
[`flux_balance`](@ref).
"""
function optimized_constraints(
    constraints::C.ConstraintTreeElem;
    modifications = [],
    output = constraints,
    kwargs...,
)
    om = optimization_model(constraints; kwargs...)
    for m in modifications
        m(om)
    end
    J.optimize!(om)
    is_solved(om) ? C.constraint_values(output, J.value.(om[:x])) : nothing
end

export optimized_constraints
