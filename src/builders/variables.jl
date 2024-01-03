
"""
$(TYPEDSIGNATURES)

Allocate a variable for each item in a constraint tree (or any other kind of
tree) and return a tree with variables bounded by the `makebound` function
which converts a given value into a bound for the corresponding variable.
"""
variables_for(makebound, ts::C.Tree) =
    let var_idx = 0
        C.tree_map(ts, C.Constraint) do x
            var_idx += 1
            C.variable(idx = var_idx, bound = makebound(x))
        end
    end
