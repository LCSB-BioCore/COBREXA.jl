
sum_objective(x) = C.Constraint(sum(c.value for c=x))
sum_objective(x::ConstraintTree) = sub_objective(values(x))

squared_error_objective(x) = C.Constraint(sum(c.value * c.value for c=x, init=zero(C.Value)))
squared_error_objective(x::ConstraintTree) = squared_error_objective(values(x))

squared_error_objective(constraints::Vector,target::Vector) =
    C.Constraint(
        sum(let tmp=(c.value - t); tmp*tmp; end for (c,t)=zip(constraints, target)))

squared_error_objective(constraints::ConstraintTree, target) =
    C.Constraint(
        sum(let tmp=(c.value - target[k]); tmp*tmp; end for (k,c)=C.elems(constraints) if haskey(target, k))
    )
