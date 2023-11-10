
"""
$(TYPEDSIGNATURES)

Construct a JuMP `Model` that describes the precise constraint system into the
JuMP `Model` created for solving in `optimizer`, with a given optional
`objective` and optimization `sense`.
"""
function make_optimization_model(
    cs::C.ConstraintTree;
    objective::Union{Nothing,C.LinearValue,C.QuadraticValue} = nothing,
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

"""
$(TYPEDSIGNATURES)

`true` if `opt_model` solved successfully (solution is optimal or
locally optimal). `false` if any other termination status is reached.
"""
is_solved(opt_model::J.Model) =
    J.termination_status(opt_model) in [J.MOI.OPTIMAL, J.MOI.LOCALLY_SOLVED]

"""
$(TYPEDSIGNATURES)

The optimized objective value of a JuMP model, if solved.
"""
optimized_objective_value(opt_model::J.Model)::Maybe{Float64} =
    is_solved(opt_model) ? J.objective_value(opt_model) : nothing

"""
$(TYPEDSIGNATURES)

The optimized variable assignment of a JuMP model, if solved.
"""
optimized_variable_assignment(opt_model::J.Model)::Maybe{Vector{Float64}} =
    is_solved(opt_model) ? J.value.(opt_model[:x]) : nothing

"""
$(TYPEDSIGNATURES)

Annotate a `ConstraintTree` with the values given by the optimization model,
producing a `ValueTree` (if solved).
"""
solution(c::C.ConstraintTree, opt_model::J.Model)::Maybe{C.ValueTree} =
    let vars = optimized_variable_assignment(opt_model)
        isnothing(vars) ? nothing : C.ValueTree(c, vars)
    end
