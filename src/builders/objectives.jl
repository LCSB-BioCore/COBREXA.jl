
sum_objectve(x) = C.Constraint(sum(C.value.(x), init = zero(C.LinearValue)))

sum_objective(x::C.ConstraintTree) = squared_error_objective(values(x))

squared_sum_objective(x) =
    C.Constraint(sum(C.squared.(C.value.(x)), init = zero(C.QuadraticValue)))

squared_sum_objective(x::C.ConstraintTree) = squared_sum_objective(values(x))

squared_error_objective(constraints::C.ConstraintTree, target) = C.Constraint(
    sum(let tmp = (C.value(c) - target[k])
            tmp * tmp
        end for (k, c) in constraints if haskey(target, k)),
)

# TODO use `mergewith` to do this reasonably (add it to ConstraintTrees)
