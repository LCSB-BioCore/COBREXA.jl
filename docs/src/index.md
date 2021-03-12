```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# Constraint-Based Reconstruction and EXascale Analysis

# Functions

A full reference to all functions is given here:

```@contents
Pages = ["functions.md"]
```

# How to contribute?

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

## Installation

To install this package: `] add https://github.com/stelmo/CobraTools.jl`.

Some of the optional features used in this package require external programs and/or data to be available. These are described below:

* The Equilibrator interface requires that the Equilibrator-API has been installed and can be accessed through Julia's PyCall package. Refer to the [Equilibrator-API website](https://gitlab.com/equilibrator/equilibrator-api) for installation instructions. Within Julia, if you can call `pyimport("equilibrator_api")` successfully, then you will be able to use the functions exposed here. To actually use the functions insert `using PyCall` in your main level script (before or after `using CobraTools`).
* To extract turnover numbers, Km, Kcat/Km and Ki from the Brenda database, you will need to download the database as a txt file [available here](https://www.brenda-enzymes.org/download_brenda_without_registration.php) (~250 MB).

The optimization solvers are implemented through `JuMP` and thus this package should be solver agnostic. All tests are conducted using `Tulip.jl` and `OSQP.jl`, but other solvers should also work (I mostly use `Gurobi.jl`). 

## Quick Example
Let's perform flux balance analysis on a constraint based model.
```@setup intro
model_location = joinpath("..","..", "models", "e_coli_core.json")
```
```@example intro
using CobraTools
using JuMP
using Tulip

model = read_model(model_location)

biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
optimizer = Tulip.Optimizer
sol = fba(model, biomass, optimizer)
```



