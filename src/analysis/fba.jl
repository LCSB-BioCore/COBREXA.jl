"""
    flux_balance_analysis(model::M, optimizer) where {M<:MetabolicModel}

Flux balance analysis solves the following problem for the input `model`:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```

Returns a solved model from [`optimize_model`](@ref).
"""
flux_balance_analysis(model::M, optimizer) where {M<:MetabolicModel} =
    optimize_model(model, optimizer; sense = MOI.MAX_SENSE)

"""
    flux_balance_analysis_vec(args...)::Union{Vector{Float64},Nothing}

A variant of FBA that returns a vector of fluxes in the same order as reactions
of the model, if the solution is found.
Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_vec(args...)::Union{Vector{Float64},Nothing}
    optmodel = flux_balance_analysis(args...)
    vars = optmodel[:x]

    termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    value.(vars)
end

"""
    flux_balance_analysis_dict(model::M, args...)::Union{Dict{String, Float64},Nothing} where {M <: MetabolicModel}

A variant of FBA that returns a dictionary assigning fluxes to reactions, if
the solution is found. Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_dict(
    model::M,
    args...,
)::Union{Dict{String,Float64},Nothing} where {M<:MetabolicModel}
    v = flux_balance_analysis_vec(model, args...)
    isnothing(v) && return nothing
    Dict(zip(reactions(model), v))
end

"""
    flux_balance_analysis(model::CobraModel, optimizer; objective_func::Union{Reaction, Array{Reaction, 1}}=Reaction[], weights=Float64[], solver_attributes=Dict{Any, Any}(), constraints=Dict{String, Tuple{Float64,Float64}}())

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
function flux_balance_analysis(
    model::CobraModel,
    optimizer;
    modifications::Union{Vector, Function} = [(model, opt_model)->nothing]
)
    # get core optimization problem
    cbm = make_optimization_model(model, optimizer)

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

    if status
        return map_fluxes(v, model)
    else
        @warn "Optimization issues occurred."
        return Dict{String,Float64}()
    end
end
