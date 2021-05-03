<div align="center">
    <img src="docs/src/assets/header.svg?maxAge=0" width="80%">
</div>

# COnstraint-Based Reconstruction and EXascale Analysis

[docs-img]:https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: http://lcsb-biocore.github.io/COBREXA.jl

[ci-img]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml/badge.svg?branch=master
[ci-url]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl/branch/master/graph/badge.svg?token=H3WSWOBD7L
[cov-url]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl

[contrib-img]: https://img.shields.io/badge/contributions-start%20here-green
[contrib-url]: https://github.com/LCSB-BioCore/COBREXA.jl/blob/master/.github/CONTRIBUTING.md

| **Documentation** | **Tests** | **Coverage** | **How to contribute?** |
|:--------------:|:-------:|:---------:|:---------:|
| [![docs-img]][docs-url] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] | [![contrib][contrib-img]][contrib-url] |

This is package provides constraint-based reconstruction and analysis tools for exa-scale metabolic models in Julia.

## How to get started

### Prerequisites and requirements

- **Operating system**: Use Linux (Debian, Ubuntu or centOS), MacOS, or Windows 10 as your operating system. `COBREXA` has been tested on these systems.
- **Julia language**: In order to use `COBREXA`, you need to install Julia 1.0 or higher. Download and follow the installation instructions for Julia [here](https://julialang.org/downloads/).
- **Hardware requirements**: `COBREXA` runs on any hardware that can run Julia, and can easily use resources from multiple computers interconnected on a network. For processing large datasets, you are required to ensure that the total amount of available RAM on all involved computers is larger than the data size.
- **Optimization solvers**: `COBREXA` uses [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) to formulate optimization problems and is compatible with all [`JuMP` supported solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers). However, to perform analysis at least one of these solvers needs to be installed on your machine. For a pure Julia implementation, we recommend [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl), but any other solver would also work.

:bulb: If you are new to Julia, it is advisable to [familiarize yourself with
the environment
first](https://docs.julialang.org/en/v1/manual/getting-started/).  Use the Julia [documentation](https://docs.julialang.org) to solve various
language-related issues, and the [Julia package manager
docs](https://julialang.github.io/Pkg.jl/v1/getting-started/) to solve
installation-related difficulties. Of course, [the Julia channel](https://discourse.julialang.org/) is another fast and easy way to find
answers to Julia specific questions.


### Installation

Using the Julia package manager to install `COBREXA` is straightforward -- after starting Julia, type:

```julia
] add COBREXA
```

> All these commands should be run from Julia at the `julia>` prompt.

Then you can load the `COBREXA` package and start using it through:

```julia
using COBREXA
```

When using `COBREXA` for the first time it may take several minutes to load, due to pre-compilation of the source code and dependencies, especially on a fresh Julia installation.

### Test the installation

If you run a non-standard platform (e.g. a customized operating systems), or if you added any modifications to the `COBREXA` source code, you may want to run the test suite to ensure that everything works as expected:

```julia
] test COBREXA
```

## Quick start guide

In addition to `COBREXA`, you also need to include a Julia package which provides an appropriate solver. One such solver is `Tulip`, which is provided by the [Tulip.jl](https://github.com/ds4dm/Tulip.jl) package.

```julia
] add Tulip
```

With the package installed and tested, let's perform simple flux balance analysis on a constraint based model.

```julia
using COBREXA
using Tulip

if !isfile("e_coli_core.xml")
  download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")
end

model = load_model("e_coli_core.xml")

sol = flux_balance_analysis_dict(model, Tulip.Optimizer)

sol["BIOMASS_Ecoli_core_w_GAM"] # 0.87
```

More functionality is described in the documentation, e.g. model construction and exa-scale analysis in pure Julia.

## Acknowledgements

`COBREXA.jl` is developed at the Luxembourg Centre for Systems Biomedicine of
the University of Luxembourg ([uni.lu/lcsb](https://www.uni.lu/lcsb)),
cooperating with the Institute for Quantitative and Theoretical Biology at the Heinrich
Heine University in Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/)).

The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://www.permedcoe.eu/)) agreement no.
951773.

<img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px">   <img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px">   <img src="docs/src/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px">
