
# Basic analysis of constraint-based models


`COBREXA.jl` supports several common analysis methods that are often used for
exploring the biological models. The currently supported methods include

- Flux balance analysis (FBA), in function [`flux_balance_analysis`](@ref)
- Flux variability analysis (FVA), in [`flux_variability_analysis`](@ref)
- Flux sampling by hit-and-run algorithm, in [`hit_and_run`](@ref)
- Parsimonious flux balance analysis (pFBA), in
  [`parsimonious_flux_balance_analysis`](@ref)

Other analysis methods are either in development and testing, or may be
specified or customized by the user. Implementing new analyses is generally
feasible; you may want to watch [the `COBREXA.jl`
repository](https://github.com/LCSB-BioCore/COBREXA.jl) for newly incoming
code, or send a feature request.

`COBREXA.jl` exports several helper functions that may help you in running
custom analysis methods:

- there is functionality to convert all types of [`MetabolicModel`](@ref)s to
  `JuMP.jl` models, mainly the function [`make_optimization_model`](@ref),
  which you may explore and analyze independently of `COBREXA.jl` using the
  tools provided by `JuMP.jl`
- there is a system of analysis modifications that allows you to easily specify
  various adjustment to the existing analysis methods

!!! tip "Notebook available!"
    Examples of running the balance functions are [available
    here](../notebooks/2_basic_usage.md).

## Optimization problem solvers

For most solving most analysis tasks, you need an optimization problem solver
that is compatible with `JuMP.jl`. You may refer to the [official list of
supported
solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers), but
we generally recommend to use either of these:

- `Tulip` (pure Julia implementation) for linear problems
- `GLPK` (based on a C library) for linear and mixed-integer problems
- `Gurobi` (based on an external library, but requires a license that is free
  for academic use) for linear, mixed-integer and quadratic problems
- `OSQP` as a free alternative to `Gurobi` for solving quadratic problems

All solvers can be installed using the Julia package manger.

## Running a single Flux balance analysis

The above methods generally accept 2 arguments: the model, and the optimizer.

In particular, having installed e.g. GLPK from the above solvers, you can run
the FBA on [your favorite *E. Coli* core
model](http://bigg.ucsd.edu/models/e_coli_core) as follows:

```
using COBREXA
m = load_model(CoreModel, "e_coli_core.xml")

using GLPK
opt_model = flux_balance_analysis(m, GLPK.Optimizer)
```

After a short while (the solver machinery usually needs to get precompiled
before the first use), you should get the `opt_model`, which is now the
optimized `JuMP.jl` model. It may print an information like this one:

```
A JuMP Model
Maximization problem with:
Variables: 95
Objective function type: JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
`JuMP.GenericAffExpr{Float64,JuMP.VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 73 constraints
`JuMP.GenericAffExpr{Float64,JuMP.VariableRef}`-in-`MathOptInterface.LessThan{Float64}`: 192 constraints
Model mode: AUTOMATIC
CachingOptimizer state: ATTACHED_OPTIMIZER
Solver name: GLPK
Names registered in the model: lbs, mb, ubs, x
```

From that, you can extract the required information with the JuMP interface,
loaded with `using JuMP`. With that,

- `objective_value(opt_model)` prints roughly `0.87`,
- `value.(opt_model[:x])` prints the vector of individual reaction fluxes.

For convenience, you can get the results nicely formatted without manually
getting them out of the optimized models:

- [`flux_balance_analysis_vec`](@ref) works like
  [`flux_balance_analysis`](@ref), but returns the vector of fluxes directly
  (in the same order as in `reactions(m)`)
- [`flux_balance_analysis_dict`](@ref) returns a dictionary with the fluxes,
  keyed by reaction identifier

## Flux variability analysis

FVA is implemented in [`flux_variability_analysis`](@ref), which returns the
usual matrix of minimal and maximal feasible flux for each reaction in the
model.

The result of calling `flux_variability_analysis(m, GLPK.Optimizer)` may look
like this (possibly showing minor numeric errors in the GLPK optimizer):

```
95×2 Array{Float64,2}:
   0.0            0.0
   6.00725        6.00725
   ⋮            
   3.64414e-13    3.17348e-13
   3.2149         3.2149
```

You can relax the optimality requirement of the reactions by specifying a wider
objective bound, getting a wider range of reaction fluxes, e.g. using
[`gamma_bounds`](@ref) (for COBRA-like γ-bound) and [`objective_bounds`](@ref)
(for a multiplicative bound around the original optimal objective).

As a result, `flux_variability_analysis(m, GLPK.Optimizer; bounds=gamma_bounds(0.8))`
will return much less constrained system:

```
95×2 Array{Float64,2}:
   0.0            0.0
   0.754299      10.1285
   ⋮            
  -4.42865        0.0
   2.57192        3.2149
```

You may additionally restrict the analysis to the list of reactions (passing
them as the second argument), and request a dictionary of the resulting fluxes,
with [`flux_variability_analysis_dict`](@ref).

## Parsimonous flux balance analysis

pFBA requires a solver that can handle quadratic problems. You may use `OSQP`
or `Gurobi`.

Otherwise, the function behaves as `flux_balance_analysis`:

- `parsimonous_flux_balance_analysis(m, OSQP.Optimizer)` will get you a
  `JuMP.jl` model optimized to a slightly more realistic (parsimonous) optimum
  than plain FBA,
- [`parsimonous_flux_balance_analysis_vec`](@ref) will return the fluxes in a
  vector,
- [`parsimonous_flux_balance_analysis_dict`](@ref) will return a reaction-keyed
  dictionary.

## Flux sampling

For the [`hit_and_run`](@ref), you need a previously optimized and constrained
model from another analysis function, such as [`flux_balance_analysis`](@ref),
or created by [`make_optimization_model`](@ref). You may need to carefully
choose the number of iterations and sample sizes to match your model; see the
documentation of [`hit_and_run`](@ref) for details.

Other than that, you can run the sampling for 100 thousand iterations as such:
```
hit_and_run(100_000, make_optimization_model(m, GLPK.Optimizer))
```

You should receive a matching flux sample with the (default) 1000 samples in a
matrix that may look like this one:
```
95×1000 Array{Float64,2}:
   0.0           0.0         …   0.0
   7.82669       9.38895         3.30653
   7.13016       4.36813         9.64434
  -0.290925     -9.3037         -0.0908829
  24.1294       17.4794          0.0511032
   ⋮                         ⋱  
 -16.243       -37.4763         -5.57301
   0.0           0.0             0.0
  -0.310819     -1.20057e-7     -2.13126
   5.71597e-5    0.00990677      0.692399
```
