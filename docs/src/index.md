# CobraTools.jl
*CobraTools is a Julia package for constraint based reconstruction and analysis of metabolic models.*

## Contents
```@contents
    Pages = [
        "io.md",
        "model_construction.md",
        "basic_analysis.md",
        "sampling_tools.md",
        "thermo_tools.md",
        "brenda_tools.jl"
    ]
    Depth=1
```

## Installation

To install this package: `] add CobraTools`.

Some of the features used in this package require external programs to be installed. These are described below:

* To import Matlab formatted models you will need `MATLAB.jl` installed and a working Matlab. Refer to the documentation of `MATLAB.jl` for installation instructions.
* The Equilibrator interface requires that the Equilibrator-API has been installed and can be accessed through Julia's PyCall package. Refer to the [Equilibrator-API website](https://gitlab.com/equilibrator/equilibrator-api) for installation instructions. Within Julia, if you can call `pyimport("equilibrator_api")` successfully, then you will be able to use the functions exposed here.
* To extract turnover numbers, Km, Kcat/Km and Ki from the Brenda database, you will need to download the database as a txt file [available here](https://www.brenda-enzymes.org/download_brenda_without_registration.php) (~250 MB).

The optimization solvers are implemented through `JuMP` and thus this package should be solver agnostic. All tests are conducted using `Gurobi.jl` but other solvers should work. 

## Quick Example
```@example
a = 1
b = 2
a+b
```

