"""
$(TYPEDSIGNATURES)

Run parsimonious flux balance analysis (pFBA) on the `model`. In short, pFBA
runs two consecutive optimization problems. The first is traditional FBA:
```
max cᵀx = μ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
And the second is a quadratic optimization problem:
```
min Σᵢ xᵢ²
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     μ = μ⁰
```
Where the optimal solution of the FBA problem, μ⁰, has been added as an
additional constraint. See "Lewis, Nathan E, Hixson, Kim K, Conrad, Tom M,
Lerman, Joshua A, Charusanti, Pep, Polpitiya, Ashoka D, Adkins, Joshua N,
Schramm, Gunnar, Purvine, Samuel O, Lopez-Ferrer, Daniel, Weitz, Karl K, Eils,
Roland, König, Rainer, Smith, Richard D, Palsson, Bernhard Ø, (2010) Omic data
from evolved E. coli are consistent with computed optimal growth from
genome-scale models. Molecular Systems Biology, 6. 390. doi:
accession:10.1038/msb.2010.47" for more details.

pFBA gets the model optimum by standard FBA (using
[`flux_balance_analysis`](@ref) with `optimizer` and `modifications`), then
finds a minimal total flux through the model that still satisfies the (slightly
relaxed) optimum. This is done using a quadratic problem optimizer. If the
original optimizer does not support quadratic optimization, it can be changed
using the callback in `qp_modifications`, which are applied after the FBA. See
the documentation of [`flux_balance_analysis`](@ref) for usage examples of
modifications.

Thhe optimum relaxation sequence can be specified in `relax` parameter, it
defaults to multiplicative range of `[1.0, 0.999999, ..., 0.99]` of the original
bound.

Returns an optimized model that contains the pFBA solution (or an unsolved model
if something went wrong).

# Example
```
model = load_model("e_coli_core.json")
optmodel = parsimonious_flux_balance_analysis(model, biomass, Gurobi.Optimizer)
value.(solution[:x])  # extract the flux from the optimizer
```
"""
function parsimonious_flux_balance_analysis(
    model::AbstractMetabolicModel,
    optimizer;
    modifications = [],
    qp_modifications = [],
    relax_bounds = [1.0, 0.999999, 0.99999, 0.9999, 0.999, 0.99],
)
    # Run FBA
    opt_model = flux_balance_analysis(model, optimizer; modifications = modifications)
    is_solved(opt_model) || return opt_model # FBA failed

    # get the objective
    Z = objective_value(opt_model)
    original_objective = objective_function(opt_model)

    # prepare the model for pFBA
    for mod in qp_modifications
        mod(model, opt_model)
    end

    # add the minimization constraint for total flux
    v = opt_model[:x] # fluxes
    @objective(opt_model, Min, sum(dot(v, v)))

    for rb in relax_bounds
        lb, ub = objective_bounds(rb)(Z)
        @models_log @info "pFBA step relaxed to [$lb,$ub]"
        @constraint(opt_model, pfba_constraint, lb <= original_objective <= ub)

        optimize!(opt_model)
        is_solved(opt_model) && break

        delete(opt_model, pfba_constraint)
        unregister(opt_model, :pfba_constraint)
    end

    return opt_model
end

"""
$(TYPEDSIGNATURES)

Perform parsimonious flux balance analysis on `model` using `optimizer`.
Returns a vector of fluxes in the same order as the reactions in `model`.
Arguments are forwarded to [`parsimonious_flux_balance_analysis`](@ref)
internally.

This function is kept for backwards compatibility, use [`flux_vector`](@ref)
instead.
"""
parsimonious_flux_balance_analysis_vec(model::AbstractMetabolicModel, args...; kwargs...) =
    flux_vector(model, parsimonious_flux_balance_analysis(model, args...; kwargs...))

"""
$(TYPEDSIGNATURES)

Perform parsimonious flux balance analysis on `model` using `optimizer`.
Returns a dictionary mapping the reaction IDs to fluxes. Arguments are
forwarded to [`parsimonious_flux_balance_analysis`](@ref) internally.

This function is kept for backwards compatibility, use [`flux_dict`](@ref)
instead.
"""
parsimonious_flux_balance_analysis_dict(model::AbstractMetabolicModel, args...; kwargs...) =
    flux_dict(model, parsimonious_flux_balance_analysis(model, args...; kwargs...))
