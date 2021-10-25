"""
    minimize_metabolic_adjustment(
        model::MetabolicModel,
        flux_ref::Vector{Float64},
        optimizer;
        modifications = [],
    )

Run minimization of metabolic adjustment (MOMA) on `model` with respect to `flux_ref`, which
is a vector of fluxes in the order of `reactions(model)`. MOMA find the shortest Euclidian
distance between `flux_ref` and `model` with `modifications`:
```
min Σᵢ (xᵢ - flux_refᵢ)²
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
Note, the problem has a quadratic constraint, so a QP solver is required. See "Daniel,
Vitkup & Church, Analysis of Optimality in Natural and Perturbed Metabolic Networks,
Proceedings of the National Academy of Sciences, 2002" for more details.

Returns an optimized model that contains the resultant nearest flux.

# Example
```
model = load_model("e_coli_core.json")
flux_ref = flux_balance_analysis_vec(model, Gurobi.Optimizer)
optmodel = minimize_metabolic_adjustment(
    model,
    flux_ref,
    Gurobi.Optimizer;
    modifications = [knockout("PFL")], # find flux of mutant that is closest to the wild type (reference) model
    )
value.(solution[:x])  # extract the flux from the optimizer
```
"""
function minimize_metabolic_adjustment(
    model::MetabolicModel,
    flux_ref::Vector{Float64},
    optimizer;
    modifications = [],
)
    opt_model = make_optimization_model(model, optimizer)

    # prepare the model for MOMA
    for mod in modifications
        mod(model, opt_model)
    end

    # moma objective
    v = opt_model[:x] # fluxes
    @objective(opt_model, Min, sum((v[i] - flux_ref[i])^2 for i=1:n_reactions(model)))

    optimize!(opt_model)

    is_solved(opt_model) || return nothing # MOMA failed

    return opt_model
end

"""
    minimize_metabolic_adjustment(
        model::MetabolicModel,
        flux_ref::Dict{String, Float64},
        optimizer;
        modifications = [],
    )

A variant of [`minimize_metabolic_adjustment`](@ref) that accepts a dictionary mapping
reaction ids to reference fluxes instead of a vector.
"""
minimize_metabolic_adjustment(
    model::MetabolicModel,
    flux_ref::Dict{String, Float64},
    optimizer;
    modifications = [],
) = minimize_metabolic_adjustment(
    model,
    [flux_ref[k] for k in reactions(model)],
    optimizer;
    modifications
)

"""
    minimize_metabolic_adjustment_vec(args...; kwargs...)

Perform minimization of metabolic adjustment (MOMA) and return a vector of fluxes in the
same order as the reactions in `model`. Arguments are forwarded to
[`minimize_metabolic_adjustment`](@ref) internally.
"""
function minimize_metabolic_adjustment_vec(args...; kwargs...)
    opt_model = minimize_metabolic_adjustment(args...; kwargs...)

    isnothing(opt_model) && return nothing

    return value.(opt_model[:x])
end

"""
    minimize_metabolic_adjustment_dict(args...; kwargs...)

Perform minimization of metabolic adjustment (MOMA) and return a dictionary mapping the
reaction IDs to fluxes. Arguments are forwarded to [`minimize_metabolic_adjustment`](@ref)
internally.
"""
function minimize_metabolic_adjustment_dict(model::MetabolicModel, args...; kwargs...)
    opt_fluxes = minimize_metabolic_adjustment_vec(model, args...; kwargs...)

    isnothing(opt_fluxes) && return nothing

    return Dict(zip(reactions(model), opt_fluxes))
end
