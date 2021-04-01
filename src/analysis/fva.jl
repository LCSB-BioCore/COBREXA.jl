"""
    flux_variability_analysis(
        model::LM,
        reactions::Vector{Int},
        optimizer,
        workers = [myid()];
        gamma::AbstractFloat = 1.0,
    )::Matrix{Float64} where {LM<:MetabolicModel}

# Flux variability analysis (FVA)

FVA solves the pair of optimization problems in `model` for each flux xᵢ listed
in `reactions`
```
min/max xᵢ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     cᵀx ≥ γ Z₀
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution to the associated
FBA problem.

Internally uses the specified JuMP-compatible `optimizer`, and the work is
scheduled in parallel on `workers`.

Returns a matrix of minima and maxima of size (length(reactions),2).
"""
function flux_variability_analysis(
    model::LM,
    reactions::Vector{Int},
    optimizer,
    workers = [myid()];
    gamma::AbstractFloat = 1.0,
)::Matrix{Float64} where {LM<:MetabolicModel}

    if any(reactions .< 1) || any(reactions .> n_reactions(model))
        throw(DomainError(reactions, "Index exceeds number of reactions."))
    end

    optimization_model = flux_balance_analysis(model, optimizer)
    Z0 = JuMP.objective_value(optimization_model)
    optimization_model = nothing # we won't need this one anymore, so free the memory

    # store a JuMP optimization model at all workers
    save_model = :(
        begin
            optmodel = COBREXA.make_optimization_model($model, $optimizer)
            x = optmodel[:x]
            COBREXA._FVA_add_constraint(optmodel, $(objective(model)), x, $Z0, $gamma)
            optmodel
        end
    )
    map(fetch, save_at.(workers, :cobrexa_parfva_model, Ref(save_model)))
    save_model = nothing # this has some volume, free it again

    # schedule FVA parts parallely using pmap
    fluxes = dpmap(
        rid -> :(COBREXA._FVA_optimize_reaction(cobrexa_parfva_model, $rid)),
        CachingPool(workers),
        [-reactions reactions],
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end

"""
    flux_variability_analysis(
        model::LM,
        optimizer;
        gamma::AbstractFloat = 1.0,
    ) where {LM<:MetabolicModel}

A simpler version of FVA that maximizes and minimizes all reactions in the model.
"""
function flux_variability_analysis(
    model::LM,
    optimizer;
    gamma::AbstractFloat = 1.0,
) where {LM<:MetabolicModel}
    n = n_reactions(model)
    return flux_variability_analysis(model, collect(1:n), optimizer; gamma = gamma)
end


"""
    _FVA_add_constraint(model, c, x, Z0, gamma)

Internal helper function for adding constraints to a model. Exists mainly
because for avoiding namespace problems on remote workers.
"""
function _FVA_add_constraint(model, c, x, Z0, gamma)
    JuMP.@constraint(model, c' * x ≥ gamma * Z0)
end

"""
    _FVA_get_opt(model, rid)

Helper for creating the optimized model on a remote worker, for avoiding
namespace problems.
"""
function _FVA_optimize_reaction(model, rid)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = JuMP.all_variables(model)[abs(rid)]

    JuMP.@objective(model, sense, var)
    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end

"""
    fva(model::CobraModel, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; optimum_bound=0.9999, weights=Float64[], solver_attributes=Dict{Any, Any}(), constraints=Dict{String, Tuple{Float64,Float64}}())

Run flux variability analysis (FVA) on the `model` with `objective_rxn(s)` and optionally specifying their `weights` (empty `weights` mean equal weighting per reaction).
It runs fba on the model once to determine the optimum of the objective.
Optionally also specify any additional flux constraints with `constraints`, a dictionary mapping reaction `id`s to tuples of (ub, lb) flux constraints.
The model is then constrained to produce objective flux bounded by `optimum_bound` from below (set to slightly less than 1.0 for stability) and each flux in the model sequentially minimized and maximized.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
The `solver_attributes` can also be specified in the form of a dictionary where each (key, value) pair will be passed to `set_optimizer_attribute(cbmodel, key, value)`.
This function builds the optimization problem from the model, and hence uses the constraints implied by the model object.
Returns two dictionaries (`fva_max` and `fva_min`) that each reaction `id`s to dictionaries of the resultant flux distributions (if solved successfully) when that `id` is optimized.

# Example
```
optimizer = Gurobi.Optimizer
atts = Dict("OutputFlag" => 0)
model = Cobraread_model("iJO1366.json")
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")
fva_max, fva_min = fva(model, biomass, optimizer; solver_attributes=atts)
```
"""
function fva(
    model::CobraModel,
    optimizer;
    optimum_bound = 0.9999,
    modifications = [(model, opt_model) -> nothing]
)
    # get core optimization problem
    cbm = make_optimization_model(model, optimizer)
    v = cbm[:x]

    # apply callbacks - user can also just put in a function
    if typeof(modifications) <: Vector
        for mod in modifications
            mod(model, cbm)
        end    
    else
        modifications(model, cbm)
    end

    optimize!(cbm)

    status = (
        termination_status(cbm) == MOI.OPTIMAL ||
        termination_status(cbm) == MOI.LOCALLY_SOLVED
    )

    fva_min = Dict{String,Dict{String,Float64}}()
    fva_max = Dict{String,Dict{String,Float64}}()

    if !status
        @warn "Error getting the initial optimum, aborting "
        return fva_max, fva_min
    end

    # Now do FVA
    λ = objective_value(cbm) # objective value
    @constraint(
        cbm,
        optimum_bound * λ <= sum(opt_weights[i] * v[i] for i in objective_indices) <= λ
    ) # for stability

    for i = 1:length(v)
        @objective(cbm, Max, v[i])
        optimize!(cbm)
        status = (
            termination_status(cbm) == MOI.OPTIMAL ||
            termination_status(cbm) == MOI.LOCALLY_SOLVED
        )
        if status
            fva_max[model.reactions[i].id] = map_fluxes(v, model)
        else
            @warn "Error maximizing index: $i with error $(termination_status(cbm))"
        end

        @objective(cbm, Min, v[i])
        optimize!(cbm)
        status = (
            termination_status(cbm) == MOI.OPTIMAL ||
            termination_status(cbm) == MOI.LOCALLY_SOLVED
        )
        if status
            fva_min[model.reactions[i].id] = map_fluxes(v, model)
        else
            @warn "Error minimizing index: $i with error $(termination_status(cbm))"
        end
    end

    return fva_max, fva_min
end
