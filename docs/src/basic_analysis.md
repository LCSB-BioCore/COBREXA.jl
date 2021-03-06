# Optimization Based Analysis
A selection of standard COBRA functions have been implemented to make basic model analysis more convenient.

## Flux balance analysis (FBA)
Flux balance analysis solves the linear program,
```math
\begin{aligned}
& \underset{v}{\text{max}}
& \sum_{i}^I {w_i}{v_i} \\
& \text{s. t.}
& Sv = 0 \\
& & v_{\text{LB}} \leq v \leq v_{\text{UB}} \\
\end{aligned}
```
using any `JuMP` compatible solver. Typically ``I`` is a singleton set that only includes the index of the biomass objective function, ``v_\mu``, and weight, ``w_\mu=1``.
```@docs
fba
```
Here, we use `Tulip.jl`, a pure Julia interior point linear program solver, with the `fba` function from `CobraTools.jl`.
```@setup fba
model_location = joinpath("..","..", "models", "e_coli_core.json")
```
```@example fba
using CobraTools
using JuMP
using Tulip

model = read_model(model_location)

biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
optimizer = Tulip.Optimizer
sol = fba(model, biomass, optimizer)
```
Note that the result of the optimization problem, `sol` maps fluxes to reaction `id`s in the model, to simplify down stream analysis.
Since the fluxes are mapped to a dictionary, it makes it simple to export them as a JSON file for visualization in, e.g. [Escher](https://escher.github.io/#/).
```@example fba
using JSON

open("fluxes.json", "w") do io
    JSON.print(io, sol)
end
```
## Parsimonious FBA
Parsimonious FBA (pFBA) solves a two stage optimization problem. First, a classic FBA problem is solved to identify the unique maximum of the objective. 
However, it should be noted that the fluxes from FBA are not unique (i.e. many fluxes may yield the objective optimum). 
To yield a unique set of fluxes, and remove internal futile cycles, a secondary quadratic problem is imposed *ad hoc*. 
Suppose that FBA has found the optimum of ``v_\mu = \mu``, pFBA then solves,
```math
\begin{aligned}
& \underset{v}{\text{max}}
& \sum_{i} {v_i^2} \\
& \text{s. t.}
& Sv = 0 \\
& & v_{\text{LB}} \leq v \leq v_{\text{UB}} \\
& v_\mu = \mu
\end{aligned}
```
again using any JuMP compatible solver(s). If multiple solvers are given, the first solver is used to solve the LP, and the second solver the QP, otherwise the same solver is used to solve both problems.
This is useful if the QP solver does not handle the LP problem well, as with OSQP. 

An alternative, related formulation of this idea exists, called "CycleFreeFlux". 
This replaces the quadratic formulation with an L1 (taxi cab) norm. While this should also remove futile cycles, it doesn't have the same uniqueness qualities and doesn't
really have much benefit beyond only solving two linear programs. See [Building your own optimization analysis script](@ref) for ideas about how to implement this yourself if you really need it.
```@docs
pfba
```
Here, we use `Tulip.jl` followed by `OSQP.jl`, with the `pfba` function from `CobraTools.jl`. Note that `OSQP.jl` has iffy performance, and is only included here because it is open source. We strongly recommend that a commercial solver, e.g. `Gubobi.jl` be used to simplify your user experience.
```@example fba
using OSQP

atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) # set solver attributes for QSQP
sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts))
```
## Flux variability analysis (FVA)

## Building your own optimization analysis script
CobraTools.jl makes it simple to access the 

```@docs
map_fluxes
```