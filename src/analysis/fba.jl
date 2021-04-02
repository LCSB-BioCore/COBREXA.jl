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
function flux_balance_analysis_vec(args...; kwargs...)::Union{Vector{Float64},Nothing}
    optmodel = flux_balance_analysis(args...; kwargs...)

    JuMP.termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    value.(optmodel[:x])
end

"""
    flux_balance_analysis_dict(model::M, args...)::Union{Dict{String, Float64},Nothing} where {M <: MetabolicModel}

A variant of FBA that returns a dictionary assigning fluxes to reactions, if
the solution is found. Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_dict(
    model::M,
    args...;
    kwargs...,
)::Union{Dict{String,Float64},Nothing} where {M<:MetabolicModel}
    v = flux_balance_analysis_vec(model, args...; kwargs...)
    isnothing(v) && return nothing
    Dict(zip(reactions(model), v))
end

"""
    flux_balance_analysis(model::StandardModel, optimizer; modifications)

Run flux balance analysis (FBA) on the `model` optionally specifying `modifications` to the problem.
These modifications can be entered as an array of modifications, or a single modification.
Leave this keyword argument out if you do not want to modify the problem.
See [`modify_constraint`](@ref), [`modify_solver_attribute`](@ref),[`modify_objective`](@ref), and [`modify_sense`](@ref)
for possible modifications.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
Returns a solved JuMP model.

# Example
```
optimizer = Gurobi.Optimizer
model = CobraTools.read_model("e_coli_core.json")
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
solved_model = fba(model, optimizer; modifications=[modify_objective(biomass)])
```
"""
function flux_balance_analysis(
    model::StandardModel,
    optimizer;
    modifications = [(model, opt_model) -> nothing],
)
    # get core optimization problem
    cbm = make_optimization_model(model, optimizer)

    # apply callbacks
    if typeof(modifications) <: Vector # many modifications
        for mod in modifications
            mod(model, cbm)
        end
    else # single modification
        modifications(model, cbm)
    end

    JuMP.optimize!(cbm)

    return cbm
end
