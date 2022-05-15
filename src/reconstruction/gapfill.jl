"""
gapfill(
    model::MetabolicModel,
    universal_reactions::Vector{Reaction},
    objective_id_bounds::Tuple{String, Float64, Float64},
    optimizer;
    modifications=[], 
    weights = fill(1.0, length(universal_reactions)),
)

Return the reaction ids of `universal_reactions` that should be added to `model`
so that the model can carry flux through its objective function. The objective
is specified through `objective_id_bounds`, which is a tuple of `(id, lb, ub)`.
If a feasible solution can be found, then the objective flux must be between
`lb` and `ub`. Optionally, specify `weights` that can be used to bias the
reactions found through solving the underlying mixed integer program. This gap
filling algorithm is based on the one introduced in *Reed, Jennifer L., et al.
"Systems approach to refining genome annotation." Proceedings of the National
Academy of Sciences (2006)*. Briefly, the algorithm is:
```
min     ∑ wᵢ * yᵢ
s.t.    S * x = 0
        xₗ ≤ x ≤ xᵤ ∀ model reactions 
        y * xₗ ≤ x ≤ y * xᵤ ∀ universal reactions 
        lb ≤ x_obj ≤ ub
        y ∈ {0, 1}
```
"""
function gapfill(
    model::MetabolicModel,
    universal_reactions::Vector{Reaction},
    objective_id_bounds::Tuple{String, Float64, Float64},
    optimizer;
    modifications=[], 
    weights = fill(1.0, length(universal_reactions)),
)
    # constraints from model to be gap filled
    S_model = stoichiometry(model)
    mids = metabolites(model)
    lbs, ubs = bounds(model)

    # set objective bounds
    obj_idx = first(indexin([objective_id_bounds[1]], reactions(model)))
    lbs[obj_idx] = objective_id_bounds[2]
    ubs[obj_idx] = objective_id_bounds[3]

    # constraints from universal reactions that can fill gaps
    n_universal_reactions = length(universal_reactions)
    n_rxns_model = size(S_model, 2)
    S_universal, lbs_universal, ubs_universal = COBREXA._universal_stoichiometry(universal_reactions, mids)

    # optimization problem is a MILP 
    n_vars = size(S_model, 2) + size(S_universal, 2)
    S = [
            [
                S_model
                spzeros(size(S_universal, 1) - size(S_model, 1), size(S_model, 2))
            ] S_universal
    ]
    bal = [
        balance(model)
        spzeros(size(S_universal, 1)-size(S_model, 1))
    ]

    opt_model = JuMP.Model(optimizer)
    @variable(opt_model, x[1:n_vars]) # fluxes
    @variable(opt_model, y[1:n_universal_reactions], Bin) # indicators
    @constraint(opt_model, mb, S * x .== bal) # mass balance
    # flux bounds of model
    @constraint(opt_model, lbs, lbs .<= x[1:n_rxns_model])
    @constraint(opt_model, ubs, x[1:n_rxns_model] .<= ubs)
    # flux bounds of universal reactions with indicators
    @constraint(opt_model, lbs_universal, lbs_universal .* y .<= x[(n_rxns_model+1):end]) 
    @constraint(opt_model, ubs_universal, x[(n_rxns_model+1):end] .<= ubs_universal .* y) 

    @objective(opt_model, Min, sum(weights .* y))

    for mod in modifications
        mod(opt_model, model)
    end

    optimize!(opt_model)

    [universal_reactions[idx].id for idx in findall(value.(y) .> 0)]
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
    metabolite_id_order,
)
    rows = Int[]
    cols = Int[]
    vals = Int[]
    lbs = zeros(length(universal_reactions))
    ubs = zeros(length(universal_reactions))
    met_id_order_lu = Dict(zip(metabolite_id_order, 1:length(metabolite_id_order)))
    n_midxs = length(met_id_order_lu)

    for (col, rxn) in enumerate(universal_reactions)
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

    return sparse(rows, cols, vals, n_midxs, length(universal_reactions)), lbs, ubs
end