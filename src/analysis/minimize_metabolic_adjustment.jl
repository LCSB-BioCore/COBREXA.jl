"""
$(TYPEDSIGNATURES)

Run minimization of metabolic adjustment (MOMA) on `model` with respect to
`flux_ref`, which is a vector of fluxes in the order of `reactions(model)`.
MOMA finds the shortest Euclidian distance between `flux_ref` and `model` with
`modifications`:
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

Returns an optimized model that contains the resultant nearest flux.

# Example
```
model = load_model("e_coli_core.json")
flux_ref = flux_balance_analysis_vec(model, Gurobi.Optimizer)
optmodel = minimize_metabolic_adjustment(
    model,
    flux_ref,
    Gurobi.Optimizer;
    modifications = [change_constraint("PFL"; lower_bound=0, upper_bound=0)], # find flux of mutant that is closest to the wild type (reference) model
    )
value.(solution[:x])  # extract the flux from the optimizer
```
"""
minimize_metabolic_adjustment_analysis(
    model::AbstractMetabolicModel,
    flux_ref::Union{Dict{String,Float64},Vector{Float64}},
    optimizer;
    modifications = [],
    kwargs...,
) = flux_balance_analysis(
    model,
    optimizer;
    modifications = vcat([minimize_metabolic_adjustment(flux_ref)], modifications),
    kwargs...,
)

"""
$(TYPEDSIGNATURES)

An optimization model modification that implements the MOMA in
[`minimize_metabolic_adjustment_analysis`](@ref).
"""
minimize_metabolic_adjustment(flux_ref::Vector{Float64}) =
    (model, opt_model) -> begin
        length(opt_model[:x]) == length(flux_ref) || throw(
            DomainError(
                flux_ref,
                "length of the reference flux doesn't match the one in the optimization model",
            ),
        )
        @objective(opt_model, Min, sum((opt_model[:x] .- flux_ref) .^ 2))
    end

"""
$(TYPEDSIGNATURES)

Overload of [`minimize_metabolic_adjustment`](@ref) that works with a
dictionary of fluxes.
"""
minimize_metabolic_adjustment(flux_ref_dict::Dict{String,Float64}) =
    (model, opt_model) ->
        minimize_metabolic_adjustment([flux_ref_dict[rid] for rid in reactions(model)])(
            model,
            opt_model,
        )

"""
$(TYPEDSIGNATURES)

Perform minimization of metabolic adjustment (MOMA) and return a vector of fluxes in the
same order as the reactions in `model`. Arguments are forwarded to
[`minimize_metabolic_adjustment`](@ref) internally.

This function is kept for backwards compatibility, use [`flux_vector`](@ref)
instead.
"""
minimize_metabolic_adjustment_analysis_vec(model::AbstractMetabolicModel, args...; kwargs...) =
    flux_vector(model, minimize_metabolic_adjustment_analysis(model, args...; kwargs...))

"""
$(TYPEDSIGNATURES)

Perform minimization of metabolic adjustment (MOMA) and return a dictionary mapping the
reaction IDs to fluxes. Arguments are forwarded to [`minimize_metabolic_adjustment`](@ref)
internally.

This function is kept for backwards compatibility, use [`flux_vector`](@ref)
instead.
"""
minimize_metabolic_adjustment_analysis_dict(model::AbstractMetabolicModel, args...; kwargs...) =
    flux_dict(model, minimize_metabolic_adjustment_analysis(model, args...; kwargs...))
