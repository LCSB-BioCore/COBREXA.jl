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
be run on huge computer clusters and in HPC environments. This allows for the parallization of analysis functions, reducing the
time between model construction and predictions.

In this package, you will find the usual COBRA-like functions that interface with user chosen linear programming solvers. 
We use [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface for many solvers; meaning you can plug in 
whatever compatible solver you want, such as `Tulip.jl`, `GLPK.jl`, `OSQP.jl`, and `Gurobi`, and things will `just work`.

## Quick start guide

You can install COBREXA using the Julia package manager as usual, by pressing `]` to
switch to Packaging environment, and typing:
```
add COBREXA
```
Note, you will need to have at least a linear programming solver installed as well
to use the analysis functions. For a pure Julia solution we recommend `Tulip.jl`,
but any solver supported by `JuMP` will also work, e.g. `Gurobi.jl`.

With the package installed, let's perform a simple flux balance analysis on a
constraint based model:

```@example intro
using COBREXA
using JuMP
using Tulip

if !isfile("e_coli_core.json")
  download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
end

model = read_model("e_coli_core.json")

biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
optimizer = Tulip.Optimizer
sol = fba(model, biomass, optimizer)
```

## Tutorials

**Note:** All tutorials assume that you have downloaded the `e_coli_core.json`
model as above.

```@contents
Pages = [
    "background.md",
    "tutorials.md",
    "functions.md",
    "howToContribute.md",
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
