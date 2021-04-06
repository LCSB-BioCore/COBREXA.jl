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
    Z0 = COBREXA.JuMP.objective_value(optimization_model)
    optimization_model = nothing # we won't need this one anymore, so free the memory

    # store a JuMP optimization model at all workers
    save_model = :(
        begin
            optmodel = COBREXA.make_optimization_model($model, $optimizer)
            COBREXA._FVA_add_constraint(
                optmodel,
                $(objective(model)),
                optmodel[:x],
                $Z0,
                $gamma,
            )
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
    COBREXA.JuMP.@constraint(model, c' * x ≥ gamma * Z0)
end

"""
    _FVA_get_opt(model, rid)

Helper for creating the optimized model on a remote worker, for avoiding
namespace problems.
"""
function _FVA_optimize_reaction(model, rid)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = COBREXA.JuMP.all_variables(model)[abs(rid)]

    COBREXA.JuMP.@objective(model, sense, var)
    COBREXA.JuMP.optimize!(model)
    return COBREXA.JuMP.objective_value(model)
end

"""
    fva(model::StandardModel, optimizer; optimum_bound=1.0-DEFAULT_FVA_TOL, modifications)

Run flux variability analysis (FVA) on the `model` (of type `StandardModel`). 
Optionally specifying problem modifications like in [`flux_balance_analysis`](@ref).
This algorithm runs FBA on the model to determine the optimum of the objective.
This optimum then constrains subsequent problems, where `optimum_bound` can be used to 
relax this constraint as a fraction of the FBA optimum, e.g. 
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
Note, this function only runs serially. 
Consider using a different model type for parallel implementations. 
Returns two dictionaries (`fva_max` and `fva_min`) that map each reaction `id` to dictionaries of the resultant flux distributions when that `id` is optimized.

See also: [`LinearModel`](@ref)

# Example
```
optimizer = Gurobi.Optimizer
model = Cobraread_model("iJO1366.json")
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")
fva_max, fva_min = fva(model, biomass, optimizer; solver_attributes=atts)
```
"""
function flux_variability_analysis(
    model::StandardModel,
    optimizer;
    optimum_bound = 1.0-DEFAULT_FVA_TOL,
    modifications = [(model, opt_model) -> nothing],
)
    # Run FBA
    opt_model = flux_balance_analysis(model, optimizer; modifications = modifications)

    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        (return nothing, nothing)

    fva_min = Dict{String,Union{Nothing,Dict{String,Float64}}}()
    fva_max = Dict{String,Union{Nothing,Dict{String,Float64}}}()

    # Now do FVA
    v = opt_model[:x]

    λ = COBREXA.JuMP.objective_value(opt_model) # objective value
    λmin = min(optimum_bound*λ, λ * 1.0 / optimum_bound)
    λmax = max(optimum_bound*λ, λ * 1.0 / optimum_bound)

    COBREXA.JuMP.@constraint(
        opt_model,
        λmin <=
        COBREXA.JuMP.objective_function(opt_model) <=
        λmax # in case there is a negative bound
    )

    for i = 1:length(v)
        COBREXA.JuMP.@objective(opt_model, Max, v[i])
        COBREXA.JuMP.optimize!(opt_model)
        status = (
            COBREXA.JuMP.termination_status(opt_model) == MOI.OPTIMAL ||
            COBREXA.JuMP.termination_status(opt_model) == MOI.LOCALLY_SOLVED
        )
        if status
            fva_max[model.reactions[i].id] =
                Dict(zip(reactions(model), value.(opt_model[:x]) ))
        else
            @warn "Error maximizing index: $i with error $(termination_status(opt_model))"
            fva_max[model.reactions[i].id] = nothing
        end

        @objective(opt_model, Min, v[i])
        optimize!(opt_model)
        status = (
            termination_status(opt_model) == MOI.OPTIMAL ||
            termination_status(opt_model) == MOI.LOCALLY_SOLVED
        )
        if status
            fva_min[model.reactions[i].id] =
                Dict(zip(reactions(model), value.(opt_model[:x]) ))
        else
            @warn "Error minimizing index: $i with error $(termination_status(opt_model))"
            fva_min[model.reactions[i].id] = nothing
        end
    end

    return fva_max, fva_min
end
