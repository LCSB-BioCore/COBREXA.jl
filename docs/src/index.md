```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# COnstraint-Based Reconstruction and EXascale Analysis

__CO__nstraint-__B__ased __R__econstruction and __EX__ascale __A__nalysis (COBREXA) provides the
necessary tools for working with exa-scale constraint-based metabolic models. Modern constraint-based 
metabolic models of complex organisms or communities typically involve the construction and analysis of systems 
comprised of billions of coupled reactions and metabolites. Current constraint-based reconstruction and analysis 
(COBRA) tools, e.g. Cobrapy and the CobraToolbox, were not designed for these large-scale systems. The purpose of COBREXA 
is to bridge this gap. COBREXA is designed to efficiently store exa-scale models and can import models from all the
commonly used model types (`.mat`, `.xml`, `.json`). Furthermore, analysis functions are exa-scale ready, and can easily 
be run on huge computer clusters and in HPC environments. This allows for the parallelization of analysis functions, reducing the
time between model construction and predictions.

In this package, you will find the usual COBRA-like functions that interface with user chosen linear programming solvers. 
We use [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface for many solvers; meaning you can plug in 
whatever [compatible solver](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) you want, such as [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl), [`GLPK.jl`](https://github.com/jump-dev/GLPK.jl), [`OSQP.jl`](https://github.com/oxfordcontrol/OSQP.jl), and [`Gurobi.jl`](https://github.com/jump-dev/Gurobi.jl); and things will `just work`.

## Tutorials

**Note:** All tutorials assume that you have downloaded the `e_coli_core.json`
model as above.

```@contents
Pages = [
    "tutorials.md",
]
Depth = 2
```

## Functions reference

```@contents
Pages = ["functions.md"]
```

## Contribution guide

If you want to contribute, please read these guidelines first:

```@contents
Pages = ["howToContribute.md"]
```
