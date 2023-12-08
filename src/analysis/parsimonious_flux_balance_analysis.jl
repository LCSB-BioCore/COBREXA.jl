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

The `relax` parameter, which defaults to `0.999`, relaxes the objective bound of
the original bound to prevent numerical issues.

Returns a [`C.ValueTree`](@ref), or `nothing` if the solution could not be
found.

# Example
```
model = load_model("e_coli_core.json")
parsimonious_flux_balance_analysis(model, Gurobi.Optimizer)
```
"""
function parsimonious_flux_balance_analysis(
    ctmodel::C.ConstraintTree,
    optimizer;
    modifications = [],
    qp_modifications = [],
    relax = 0.999,
)
    # Run FBA
    opt_model = optimization_model(ctmodel; objective = ctmodel.objective.value, optimizer)

    for mod in modifications
        mod(ctmodel, opt_model)
    end

    J.optimize!(opt_model)

    is_solved(opt_model) || return nothing

    # get the objective
    Z = J.objective_value(opt_model)

    _ctmodel = ctmodel * :pfbaobjective^squared_sum_objective(ctmodel.fluxes)
    _ctmodel.objective.bound = relax == 1.0 ? Z : objective_bounds(relax)(Z)  # TODO currently breaks

    opt_model = X.optimization_model(
        _ctmodel;
        objective = ctmodel.pfbaobjective.value,
        optimizer,
        sense = X.J.MIN_SENSE,
    )

    # prepare the model for pFBA
    for mod in qp_modifications
        mod(model, opt_model)
    end

    J.optimize!(opt_model)

    C.ValueTree(_ctmodel, J.value.(opt_model[:x]))
end

"""
$(TYPEDSIGNATURES)

Variant that takes an [`A.AbstractFBCModel`](@ref) as input. All other arguments are forwarded.
"""
function parsimonious_flux_balance_analysis(model::A.AbstractFBCModel, optimizer; kwargs...)
    ctmodel = fbc_model_constraints(model)
    parsimonious_flux_balance_analysis(ctmodel, optimizer; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Pipe-able variant of [`parsimonious_flux_balance_analysis`](@ref).
"""
parsimonious_flux_balance_analysis(optimizer; modifications = []) =
    m -> parsimonious_flux_balance_analysis(m, optimizer; modifications)

export parsimonious_flux_balance_analysis
