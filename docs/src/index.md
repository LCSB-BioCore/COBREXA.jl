```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# Constraint-Based Reconstruction and EXascale Analysis

COBREXA provides the base toolkit for working with big models and running
tremendous amount of analyses in parallel. Its main purpose is to make the
COBRA approach easy to use on problems that require the use of huge computer
clusters and HPC environment, allowing it to approach the realistic
pre-exascale capacities.

In the package, you will find the usual COBRA-like functions that inteface to
the underlying linear programming solvers. We use
[`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface for
many solvers; you can plug in whatever compatible solver you want, such as
`Tulip.jl`, `GLPK.jl`, `OSQP.jl`, and `Gurobi`.

## Quick start guide

You can install COBREXA from Julia repositories as usual, by pressing `]` to switch to Packaging environment, and typing:
```
add COBREXA
```

With the package installed, let's perform a simple flux balance analysis on a constraint based model:

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
    "model_structure.md",
    "io.md",
    "model_construction.md",
    "basic_analysis.md",
    "sampling_tools.md",
    "external_tools.md",
    "thermodynamics.md"
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
