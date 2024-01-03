
variables_for(makebound, cs::C.ConstraintTree) =
    let var_idx = 0
        C.tree_map(cs, C.Constraint) do c
            var_idx += 1
            C.variable(idx = var_idx, bound = makebound(c))
        end
    end
