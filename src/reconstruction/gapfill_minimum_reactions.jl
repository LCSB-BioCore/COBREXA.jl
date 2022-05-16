"""
    gapfill_minimum_reactions(
        model::MetabolicModel,
        universal_reactions::Vector{Reaction},
        objective_lower_bound::Float64,
        optimizer;
        modifications=[], 
        weights = fill(1.0, length(universal_reactions)),
        objective_upper_bound = COBREXA._constants.default_reaction_bound,
        ignore_reactions = [],
        max_gaps_fillable = 1000_000,
    )
    
Return the indices of reactions in `universal_reactions` that should be added to
`model` so that the model can carry flux through its objective function, which
is bounded by `objective_lower_bound`. Optionally, specify `weights` that can be
used to bias the reactions found through solving the underlying mixed integer
program (MILP). Also, some reactions in `universal_reactions` can be ignored by
specifying their ids in `ignore_reactions`, this is useful to, e.g., restrict
which exchanges can be added. Finally, the limit the search space, it is
possible to specify the maximum number of gaps that can be filled through
`max_gaps_fillable`.

This gap filling algorithm is based on the one introduced in *Reed, Jennifer L.,
et al. "Systems approach to refining genome annotation." Proceedings of the
National Academy of Sciences (2006)*. Briefly, the algorithm find the smallest
number of reactions to add by solving the MILP:
```
min     ∑ wᵢ * yᵢ
s.t.    S * x = 0
        xₗ ≤ x ≤ xᵤ ∀ model reactions 
        y * xₗ ≤ x ≤ y * xᵤ ∀ universal reactions 
        lb ≤ objective(x) ≤ ub
        y ∈ {0, 1}
```
where `w` is the set of optional `weights`, `x` the fluxes, and `y` the indicator 
variables.
"""
function gapfill_minimum_reactions(
    model::MetabolicModel,
    universal_reactions::Vector{Reaction},
    objective_lower_bound::Float64,
    optimizer;
    modifications = [],
    weights = fill(1.0, length(universal_reactions)),
    objective_upper_bound = COBREXA._constants.default_reaction_bound,
    ignore_reactions = [],
    max_gaps_fillable = COBREXA._constants.max_gaps_fillable,
)
    # constraints from model to be gap filled
    S_model = stoichiometry(model)
    metabolite_id_order = metabolites(model)

    # constraints from universal reactions that can fill gaps
    n_universal_reactions = length(universal_reactions)
    S_universal, lbs_universal, ubs_universal = COBREXA._universal_stoichiometry(
        universal_reactions,
        metabolite_id_order;
        ignore_reactions,
    )

    # adjust the model stoichiometric matrix to account for additional metabolites if necessary
    S = [
        S_model
        spzeros(size(S_universal, 1) - size(S_model, 1), size(S_model, 2))
    ]
    # adjust the balance to account for additional metabolites
    bal = [
        balance(model)
        spzeros(size(S_universal, 1) - size(S_model, 1))
    ]

    #=
    First build standard flux balance type optimization problem, then add
    specific details of the gap filling algorithm, e.g. indicator constraints, etc.
    =#
    opt_model = make_optimization_model(model, optimizer; sense = COBREXA.MIN_SENSE)
    delete(opt_model, opt_model[:mb]) #  need to remove mass balances
    unregister(opt_model, :mb) # will re-use symbol

    @variable(opt_model, z[1:n_universal_reactions]) # fluxes from universal reactions 
    @variable(opt_model, y[1:n_universal_reactions], Bin) # indicators

    # objective bounds 
    @constraint(
        opt_model,
        obj_bounds,
        objective_lower_bound <= objective(model)' * opt_model[:x] <= objective_upper_bound
    )

    # flux bounds of universal reactions with indicators
    @constraint(opt_model, lbs_universal, lbs_universal .* y .<= z)
    @constraint(opt_model, ubs_universal, z .<= ubs_universal .* y)

    # combined mass balances 
    @constraint(opt_model, mb, S * opt_model[:x] + S_universal * z .== bal) # mass balance of all reactions

    # constrain the maximum number of gaps that can be filled 
    @constrain(opt_model, max_gaps, sum(y) <= max_gaps_fillable)
    
    # make new objective
    @objective(opt_model, Min, sum(weights .* y))

    for mod in modifications
        mod(opt_model, model)
    end

    optimize!(opt_model)

    findall(value.(y) .> 0)
end

"""
    _universal_stoichiometry(
        universal_reactions::Vector{Reaction},
        metabolite_id_order,
    )

A helper function that constructs the stoichiometric matrix of a set of
`universal_reactions`. The order of the metabolites is determined with
`metabolite_id_order`, so that this stoichiometric matrix can be combined with
another one.
"""
function _universal_stoichiometry(
    universal_reactions::Vector{Reaction},
    metabolite_id_order;
    ignore_reactions = [],
)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    lbs = zeros(length(universal_reactions))
    ubs = zeros(length(universal_reactions))
    met_id_order_lu = Dict(zip(metabolite_id_order, 1:length(metabolite_id_order)))
    n_midxs = length(met_id_order_lu) # account for metabolites already in model

    n_cols = 0 # counter for filtered reactions
    for (col, rxn) in
        enumerate(filter(x -> !in(x.id, ignore_reactions), universal_reactions))
        n_cols += 1
        for (mid, stoich) in rxn.metabolites
            if !haskey(met_id_order_lu, mid)
                n_midxs += 1
                met_id_order_lu[mid] = n_midxs
            end
            push!(rows, met_id_order_lu[mid])
            push!(cols, col)
            push!(vals, stoich)
        end
        lbs[col] = rxn.lb
        ubs[col] = rxn.ub
    end

    return sparse(rows, cols, vals, n_midxs, n_cols), lbs, ubs
end
