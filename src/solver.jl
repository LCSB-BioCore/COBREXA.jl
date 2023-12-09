
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
    sense = J.MAX_SENSE,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(cs)])
    isnothing(objective) || J.@objective(model, sense, C.substitute(objective, x))

    # constraints
    function add_constraint(c::C.Constraint)
        if c.bound isa Float64
            J.@constraint(model, C.substitute(c.value, x) == c.bound)
        elseif c.bound isa C.IntervalBound
            val = C.substitute(c.value, x)
            isinf(c.bound[1]) || J.@constraint(model, val >= c.bound[1])
            isinf(c.bound[2]) || J.@constraint(model, val <= c.bound[2])
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
