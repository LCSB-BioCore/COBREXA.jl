
"""
$(TYPEDSIGNATURES)

Construct a JuMP `Model` that describes the precise constraint system into the
JuMP `Model` created for solving in `optimizer`, with a given optional
`objective` and optimization `sense`.
"""
function J.Model(
    constraints::C.ConstraintTree;
    objective::Maybe{C.Value} = nothing,
    optimizer,
    sense = J.MAX_SENSE,
)
    # TODO this might better have its own name to avoid type piracy.
    model = J.Model(optimizer)
    J.@variable(model, x[1:C.var_count(cs)])

    # objectives
    if objective isa C.Value
        JuMP.@objective(model, sense, C.value_product(objective, x))
    elseif objective isa C.QValue
        JuMP.@objective(model, sense, C.qvalue_product(objective, x))
    end

    # constraints
    function add_constraint(c::C.Constraint)
        if c.bound isa Float64
            J.@constraint(model, C.value_product(c.value, x) == c.bound)
        elseif c.bound isa C.IntervalBound
            val = C.value_product(c.value, x)
            isinf(c.bound[1]) || J.@constraint(model, val >= c.bound[1])
            isinf(c.bound[2]) || J.@constraint(model, val <= c.bound[2])
        end
    end
    function add_constraint(c::C.QConstraint)
        if c.bound isa Float64
            JuMP.@constraint(model, C.qvalue_product(c.qvalue, x) == c.bound)
        elseif c.bound isa Tuple{Float64,Float64}
            val = C.qvalue_product(c.qvalue, x)
            isinf(c.bound[1]) || JuMP.@constraint(model, val >= c.bound[1])
            isinf(c.bound[2]) || JuMP.@constraint(model, val <= c.bound[2])
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

Convenience re-export of `Model` from JuMP.
"""
const Model = J.Model

"""
$(TYPEDSIGNATURES)

Convenience re-export of `optimize!` from JuMP.
"""
const optimize! = J.optimize!

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
    is_solved(opt_model) ? J.value.(model[:x]) : nothing

"""
$(TYPEDSIGNATURES)

Convenience overload for making solution trees out of JuMP models
"""
C.SolutionTree(c::C.ConstraintTree, opt_model::J.Model)::Maybe{C.SolutionTree} =
    let vars = optimized_variable_assignment(opt_model)
        isnothing(vars) ? nothing : C.SolutionTree(c, vars)
    end
