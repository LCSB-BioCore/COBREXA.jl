"""
    fluxBalanceAnalysis(model::M, optimizer) where {M<:MetabolicModel}

Flux balance analysis solves the following problem for the input `model`:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```

Returns a solved model from [`optimizeModel`](@ref).
"""
fluxBalanceAnalysis(model::M, optimizer) where {M<:MetabolicModel} =
optimizeModel(model, optimizer; sense = MOI.MAX_SENSE)

"""
    fba(model::CobraModel, optimizer; objective_func::Union{Reaction, Array{Reaction, 1}}=Reaction[], weights=Float64[], solver_attributes=Dict{Any, Any}(), constraints=Dict{String, Tuple{Float64,Float64}}())

Run flux balance analysis (FBA) on the `model` optionally specifying `objective_rxn(s)` and their `weights` (empty `weights` mean equal weighting per reaction).
Optionally also specify any additional flux constraints with `constraints`, a dictionary mapping reaction `id`s to tuples of (lb, ub) flux constraints.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work. 
The `solver_attributes` can also be specified in the form of a dictionary where each (key, value) pair will be passed to `set_optimizer_attribute(cbmodel, key, value)`.
This function builds the optimization problem from the model, and hence uses the constraints implied by the model object.
Returns a dictionary of reaction `id`s mapped to fluxes if solved successfully, otherwise an empty dictionary.

# Example
```
optimizer = Gurobi.Optimizer
atts = Dict("OutputFlag" => 0)
model = CobraTools.read_model("iJO1366.json")
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")
sol = fba(model, biomass, optimizer; solver_attributes=atts)
```
"""
function fba(
    model::CobraModel,
    optimizer;
    objective_func::Union{Reaction,Array{Reaction,1}} = Reaction[],
    weights = Float64[],
    solver_attributes = Dict{Any,Any}(),
    constraints = Dict{String,Tuple{Float64,Float64}}(),
    sense = MOI.MAX_SENSE,
)
    # get core optimization problem
    cbm, v, mb, lbcons, ubcons = makeOptimizationModel(model, optimizer, sense = sense)

    # modify core optimization problem according to user specifications
    if !isempty(solver_attributes) # set other attributes
        for (k, val) in solver_attributes
            set_optimizer_attribute(cbm, k, val)
        end
    end

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, lbcons, ubcons; lb = con[1], ub = con[2])
    end

    # if an objective function is supplied, modify the default objective
    if typeof(objective_func) == Reaction || !isempty(objective_func)
        # ensure that an array of objective indices are fed in
        if typeof(objective_func) == Reaction
            objective_indices = [model[objective_func]]
        else
            objective_indices = [model[rxn] for rxn in objective_func]
        end

        if isempty(weights)
            weights = ones(length(objective_indices))
        end
        opt_weights = zeros(length(model.reactions))

        # update the objective function tracker
        # don't update model objective function - silly thing to do
        wcounter = 1
        for i in eachindex(model.reactions)
            if i in objective_indices
                # model.reactions[i].objective_coefficient = weights[wcounter]
                opt_weights[i] = weights[wcounter]
                wcounter += 1
                # else
                # model.reactions[i].objective_coefficient = 0.0
            end
        end

        @objective(cbm, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    else # use default objective
        # automatically assigned by makeOptimizationModel
    end

    optimize!(cbm)

    status = (
        termination_status(cbm) == MOI.OPTIMAL ||
        termination_status(cbm) == MOI.LOCALLY_SOLVED
    )

    if status
        return map_fluxes(v, model)
    else
        @warn "Optimization issues occurred."
        return Dict{String,Float64}()
    end
end
