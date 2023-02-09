"""
$(TYPEDSIGNATURES)

Run minimization of metabolic adjustment (MOMA) on `model` with respect to
`reference_flux`, which is a vector of fluxes in the order of
`variables(model)`. MOMA finds the shortest Euclidian distance between
`reference_flux` and `model` with `modifications`:
```
min Σᵢ (xᵢ - flux_refᵢ)²
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
Because the problem has a quadratic objective, a QP solver is required. See
"Daniel, Vitkup & Church, Analysis of Optimality in Natural and Perturbed
Metabolic Networks, Proceedings of the National Academy of Sciences, 2002" for
more details.

Additional arguments are passed to [`flux_balance_analysis`](@ref).

Returns an optimized model that contains the feasible flux nearest to the
reference.

# Example
```
model = load_model("e_coli_core.json")
reference_flux = flux_balance_analysis_vec(model, Gurobi.Optimizer)
optmodel = minimize_metabolic_adjustment(
    model,
    reference_flux,
    Gurobi.Optimizer;
    modifications = [change_constraint("PFL"; lower_bound=0, upper_bound=0)], # find flux of mutant that is closest to the wild type (reference) model
    )
value.(solution[:x])  # extract the flux from the optimizer
```
"""
minimize_metabolic_adjustment_analysis(
    model::AbstractMetabolicModel,
    reference_flux::Union{Vector{Float64}},
    optimizer;
    kwargs...,
) = flux_balance_analysis(
    model |> COBREXA.Reconstruction.Pipes.minimize_adjustment(reference_flux),
    optimizer;
    kwargs...,
)
