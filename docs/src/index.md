```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# Constraint-Based Reconstruction and EXascale Analysis

## Installation

To install this package: `] add COBREXA`.

## Functions

A full reference to all functions is given here:

```@contents
Pages = ["functions.md"]
```

## How to contribute?

If you want to contribute, please read these guidelines first:

```@contents
Pages = ["howToContribute.md"]
```

## Contents
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
    Depth=2
```

The optimization solvers are implemented through `JuMP` and thus this package should be solver agnostic. All tests are conducted using `Tulip.jl`, `GLPK.jl`, and `OSQP.jl`, but other solvers should also work. 

## Quick Example
Let's perform flux balance analysis on a constraint based model.
```@example intro
using COBREXA
using JuMP
using Tulip

download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
model = read_model("e_coli_core.json")

biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
optimizer = Tulip.Optimizer
sol = fba(model, biomass, optimizer)
rm("e_coli_core.json") # hide
```
