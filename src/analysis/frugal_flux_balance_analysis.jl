"""
    frugal_flux_balance_analysis(
        model::MetabolicModel,
        optimizer;
        modifications = [],
        frugal_reactions = 1:n_reactions(model),
        relax_bounds = [0.999999, 0.99999, 0.9999, 0.999, 0.99],
    )

Run "frugal" flux balance analysis on the `model`. Frugal FBA is similar to pFBA
in that two sequential optimization problems are solved that remove internal
cycles. The difference being that frugal FBA replaces the quadratic (L2 norm)
objective of pFBA with an L1 norm (absolute value) objective function. The 
benefit of this is that this L1 norm problem can be converted into an LP using
standard transformations. Specifically, the first step is traditional FBA:
```
max cᵀx = μ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
And the second is an L1 norm objective optimization problem:
```
min Σᵢ |xᵢ| for i ∈ I
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     μ = μ⁰
```
Where the optimal solution of the FBA problem, μ⁰, has been added as an
additional constraint. Note that the second problem can be transformed into
and LP, hence it is more computationally efficient to perform (especially 
for large models).

By supplying `frugal_reactions` indices the reactions in set `I` can be
specified. By default, this includes all the reactions in the model. Like pBA
the objective function can be relaxed in case numerical problems occurs. The
optimum relaxation sequence can be specified in `relax` parameter, it defaults
to multiplicative range of `[1.0, 0.999999, ..., 0.99]` of the original bound.

Returns an optimized model that contains the frugal FBA solution; or `nothing`
if the optimization failed.

# Example
```
model = load_model("e_coli_core.json")
optmodel = frugal_flux_balance_analysis(model, biomass, Gurobi.Optimizer)
value.(solution[:x])  # extract the flux from the optimizer
```
"""
function frugal_flux_balance_analysis(
    model::MetabolicModel,
    optimizer;
    modifications = [],
    frugal_reactions = 1:n_reactions(model),
    relax_bounds = [1.0, 0.999999, 0.99999, 0.9999, 0.999, 0.99],
)
    if typeof(frugal_reactions) == Vector{String}
        frugal_reactions = indexin(frugal_reactions, reactions(model))
    end

    # Run FBA
    opt_model = flux_balance_analysis(model, optimizer; modifications = modifications)
    is_solved(opt_model) || return nothing # FBA failed

    # get the objective for relaxation
    Z = objective_value(opt_model)
    original_objective = COBREXA.JuMP.objective_function(opt_model)

    # add transformation variables
    t = @variable(opt_model, t[1:length(frugal_reactions)])
    @objective(opt_model, Min, sum(t))
    for (i, frugal_reaction) in enumerate(frugal_reactions)
        @constraint(opt_model, t[i] >= opt_model[:x][frugal_reaction])
        @constraint(opt_model, t[i] >= -opt_model[:x][frugal_reaction])
    end

    for rb in relax_bounds
        lb, ub = objective_bounds(rb)(Z)
        @_models_log @info "frugal FBA step relaxed to [$lb,$ub]"
        @constraint(opt_model, obj_constraint, lb <= original_objective <= ub)

        optimize!(opt_model)
        is_solved(opt_model) && break
        @warn rb
        COBREXA.JuMP.delete(opt_model, obj_constraint)
        COBREXA.JuMP.unregister(opt_model, :obj_constraint)
    end

    is_solved(opt_model) || return nothing # pFBA failed

    return opt_model
end

"""
    frugal_flux_balance_analysis_vec(args...; kwargs...)

Perform frugal flux balance analysis on `model` using `optimizer`.
Returns a vector of fluxes in the same order as the reactions in `model`.
Arguments are forwarded to [`frugal_flux_balance_analysis`](@ref)
internally.
"""
function frugal_flux_balance_analysis_vec(args...; kwargs...)
    opt_model = frugal_flux_balance_analysis(args...; kwargs...)

    isnothing(opt_model) && return nothing

    return value.(opt_model[:x])
end

"""
    frugal_flux_balance_analysis_dict(model::MetabolicModel, args...; kwargs...)

Perform frugal flux balance analysis on `model` using `optimizer`.
Returns a dictionary mapping the reaction IDs to fluxes. Arguments are
forwarded to [`frugal_flux_balance_analysis`](@ref) internally.
"""
function frugal_flux_balance_analysis_dict(model::MetabolicModel, args...; kwargs...)
    opt_fluxes = frugal_flux_balance_analysis_vec(model, args...; kwargs...)

    isnothing(opt_fluxes) && return nothing

    return Dict(zip(reactions(model), opt_fluxes))
end
