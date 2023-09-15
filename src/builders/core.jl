
import AbstractFBCModels as F
import SparseArrays: sparse

"""
$(TYPEDSIGNATURES)

Make a constraint tree that models the content of the given instance of
`AbstractFBCModel`.
"""
C.make_constraint_tree(model::F.AbstractFBCModel) = let
    rxns = Symbol.(F.reactions(model)),
    mets = Symbol.(F.metabolites(model)),
    lbs, ubs = F.bounds(model),
    stoi = F.stoichiometry(model),
    bal = F.balance(model),
    obj = F.objective(model)

    # TODO coupling
    :fluxes^C.allocate_variables(
        keys=rxns,
        bounds=zip(lbs, ubs),
    ) *
    :balance^C.make_constraint_tree(
        m => Constraint(value=Value(sparse(row)), bound=b) for
            (m, row, b) in zip(mets, eachrow(stoi), bals)
    ) *
    :objective^C.Constraint(value=C.Value(sparse(obj)))
end
