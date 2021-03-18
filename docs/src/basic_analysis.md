# Optimization Based Analysis
A selection of standard COBRA functions have been implemented to make basic model analysis more convenient.
Additionally, `COBREXA.jl` allows you to easily formulate your own optimization problems using the structure of a constraint based model.
This makes it easy to experiment with custom algorithms etc.

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
Here, we use `Tulip.jl`, a pure Julia interior point linear program solver, with the `fba` function from `COBREXA.jl`.
```@example fba
using COBREXA
using JuMP
using Tulip

model = read_model("e_coli_core.json")

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

## Solution inspection
Sometimes it is useful to investigate which reactions consume or produce a certain metabolite, or which exchange reactions are active. 
This functionality is exposed via `metabolite_fluxes` and `exchange_reactions`.
```@docs
metabolite_fluxes
exchange_reactions
```
```@example fba
consuming, producing = metabolite_fluxes(sol, model)
consuming["atp_c"]
```
## Parsimonious FBA
Parsimonious FBA (pFBA) solves a two stage optimization problem. First, a classic FBA problem is solved to identify the unique maximum of the objective. 
However, it should be noted that the fluxes from FBA are not unique (i.e. many fluxes may yield the objective optimum). 
To yield a unique set of fluxes, and remove internal futile cycles, a secondary quadratic problem is imposed *ad hoc*. 
Suppose that FBA has found the optimum of ``v_\mu = \mu``, pFBA then solves,
```math
\begin{aligned}
& \underset{v}{\text{min}}
& \sum_{i} {v_i^2} \\
& \text{s. t.}
& Sv = 0 \\
& & v_{\text{LB}} \leq v \leq v_{\text{UB}} \\
& & v_\mu = \mu
\end{aligned}
```
again using any JuMP compatible solver(s). In the `COBREXA.jl` implementation of pFBA, both the FBA and QP problem are solved internally in `pfba`, using similar input arguments as in `fba`. If multiple solvers are given, the first solver is used to solve the LP, and the second solver the QP, otherwise the same solver is used to solve both problems.
This is useful if the QP solver does not handle the LP problem well, as with OSQP. 

An alternative, related formulation of this idea exists, called "CycleFreeFlux". 
This replaces the quadratic formulation with an L1 (taxi cab) norm. While this should also remove futile cycles, it doesn't have the same uniqueness qualities and doesn't
really have much benefit beyond only solving two linear programs. See [Building your own optimization analysis script](@ref) for ideas about how to implement this yourself if you really need it.
```@docs
pfba
```
Here, we use `Tulip.jl` followed by `OSQP.jl`, with the `pfba` function from `COBREXA.jl`. Note that `OSQP.jl` has iffy performance, and is only included here because it is open source. We recommend that a commercial solver, e.g. `Gubobi.jl`, be used to simplify your user experience.
```@example fba
using OSQP

atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) # set solver attributes for QSQP
sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts))
```
## Flux variability analysis (FVA)
Flux variability analysis can also be used to investigate the degeneracy associated with flux balance analysis derived solutions (see also [Sampling Tools](@ref)).
`COBREXA.jl` exposes `fva` that sequentially maximizes and minimizes each reaction in a model subject to the constraint that each optimization problem also satisfies an initial FBA type objective optimum, below denoted by ``v_{\mu}=\mu``,
```math
\begin{aligned}
& \underset{v}{\text{max or min}}
& v_i \\
& \text{s. t.}
& Sv = 0 \\
& & v_{\text{LB}} \leq v \leq v_{\text{UB}} \\
& & v_{\mu} = \mu \\ 
\end{aligned}
```
```@docs
fva
```
## Building your own optimization analysis script
`COBREXA.jl` also makes it simple to construct customized optimization problems by making judicious use of [JuMP](https://jump.dev/). 
Convenience functions make optimization problem construction, modification and data extraction from JuMP result objects easy.
```@docs
get_core_model
build_cbm
set_bound
map_fluxes(::Array{Float64,1}, ::CobraModel)
```
```@example fba
using COBREXA
using JuMP
using Tulip

model = read_model("e_coli_core.json")
cbm, v, mb, ubs, lbs = build_cbm(model)
glucose_index = model[findfirst(model.reactions, "EX_glc__D_e")]
set_bound(glucose_index, ubs, lbs; ub=-12.0, lb=-12.0)

set_optimizer(cbm, Tulip.Optimizer)
@objective(cbm, Max, v[model[biomass]])
optimize!(cbm)    

sol = map_fluxes(v, model) 
```
