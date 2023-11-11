
squared_sum_error_objective(constraints::C.ConstraintTree, target::Dict{Symbol,Float64}) =
    C.Constraint(
        sum(
            (C.value(c) - target[k]) * (C.value(c) - target[k]) for
            (k, c) in constraints if haskey(target, k)
        ),
    )

squared_sum_objective(x::C.ConstraintTree) =
    squared_sum_error_objective(x, Dict(keys(x) .=> 0.0))

objective_bounds(tolerance) = z -> begin
    vs = (z * tolerance, z / tolerance)
    (minimum(vs), maximum(vs))
end
