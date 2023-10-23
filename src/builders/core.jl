
"""
$(TYPEDSIGNATURES)

A constraint tree that models the content of the given instance of
`AbstractFBCModel`.
"""
function metabolic_model(model::A.AbstractFBCModel)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bal = A.balance(model)
    obj = A.objective(model)

    :fluxes^C.variables(keys = rxns, bounds = zip(lbs, ubs)) *
    :flux_stoichiometry^C.ConstraintTree(
        met => C.Constraint(value = C.LinearValue(sparse(row)), bound = b) for
        (met, row, b) in zip(mets, eachrow(stoi), bal)
    ) *
    :objective^C.Constraint(C.LinearValue(sparse(obj)))
end

"""
$(TYPEDSIGNATURES)

Shortcut for allocation non-negative ("unsigned") variables. The argument
`keys` is forwarded to `ConstraintTrees.variables` as `keys`.
"""
unsigned_variables(; keys) = C.variables(; keys, bounds = Ref((0.0, Inf)))

function fluxes_in_direction(fluxes::C.ConstraintTree, direction = :forward)
    keys = Symbol[]
    for (id, flux) in fluxes
        if direction == :forward
            last(flux.bound) > 0 && push!(keys, id)
        else
            first(flux.bound) < 0 && push!(keys, id)
        end
    end
    C.variables(; keys, bounds = Ref((0.0, Inf)))
end

"""
$(TYPEDSIGNATURES)

A constraint tree that bound the values present in `signed` to be sums of pairs
of `positive` and `negative` contributions to the individual values.

Keys in the result are the same as the keys of `signed` constraints.

Typically, this can be used to create "unidirectional" fluxes
together with [`unsigned_variables`](@ref):
```
uvars = unsigned_variables(keys = collect(keys(myModel.fluxes)))

myModel = myModel +
    :fluxes_forward^uvars +
    :fluxes_reverse^uvars

myModel *=
    :direction_sums^sign_split_constraints(
        positive = myModel.fluxes_forward,
        negative = myModel.fluxes_reverse,
        signed = myModel.fluxes,
    )
```
#TODO this might go to the docs
"""
sign_split_constraints(;
    positive::C.ConstraintTree,
    negative::C.ConstraintTree,
    signed::C.ConstraintTree,
) = C.ConstraintTree(
    k => C.Constraint(
        value = s.value +
                (haskey(negative, k) ? negative[k].value : zero(typeof(s.value))) -
                (haskey(positive, k) ? positive[k].value : zero(typeof(s.value))),
        bound = 0.0,
    ) for (k, s) in signed
)
