```@raw html
<br>
<div align="center">
    <img class="docs-light-only" src="assets/header.svg?maxAge=0" width="80%">
    <img class="docs-dark-only" src="assets/header-dark.svg?maxAge=0" width="80%">
</div>
<br>
```

# Constraint-Based Reconstruction and EXascale Analysis

| **Repository** | **Tests** | **Coverage** | **How to contribute?** |
|:--------------:|:-------:|:---------:|:---------:|
| [![GitHub](https://img.shields.io/github/stars/LCSB-BioCore/COBREXA.jl?label=COBREXA.jl&style=social)](http://github.com/LCSB-BioCore/COBREXA.jl) | [![CI](https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml) | [![codecov](https://codecov.io/gh/LCSB-BioCore/COBREXA.jl/branch/master/graph/badge.svg?token=H3WSWOBD7L)](https://codecov.io/gh/LCSB-BioCore/COBREXA.jl) | [![contrib](https://img.shields.io/badge/contributions-start%20here-green)](https://github.com/LCSB-BioCore/COBREXA.jl/blob/master/.github/CONTRIBUTING.md) |


COBREXA is a toolkit for working with large constraint-based metabolic models,
and a running very large number of analysis tasks on these models in parallel.
Its main purpose is to make the methods of Constraint-based Reconstruction and
Analysis (COBRA) scale to problem sizes that require the use of huge computer
clusters and HPC environments, which allows them to be realistically applied to
pre-exascale-sized models.

In this package, you will find the usual COBRA-like functions that interface to
underlying linear programming solvers. We use
[`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface for
many solvers; you can plug in whichever compatible solver you want, including
the popular [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl),
[`GLPK.jl`](https://github.com/jump-dev/GLPK.jl),
[`OSQP.jl`](https://github.com/oxfordcontrol/OSQP.jl), and
[`Gurobi.jl`](https://github.com/jump-dev/Gurobi.jl).

```@raw html
<div align="center">
<img style="width:300px;margin:10px;border-offset:15px;border: 1px solid #eee;border-radius: 50%;padding: 10px;-webkit-border-radius: 50%;-moz-border-radius: 50%;" src="assets/output.gif" alt="history"><br>
Development history of COBREXA.jl.
</div>
```

## Quick start guide

<!--insert_quickstart-->

## Tutorials

Detailed tutorial content is [available here](tutorials.md).

```@contents
Pages = joinpath.("tutorials", filter(x -> endswith(x, ".md"), readdir("tutorials")))
Depth = 1
```

## Example notebooks and workflows

Detailed notebook content is [available here](notebooks.md).

```@contents
Pages = joinpath.("notebooks", filter(x -> endswith(x, ".md"), readdir("notebooks")))
Depth = 1
```

## Functions reference

```@contents
Pages = ["functions.md"]
```

## Contribution guide

If you wish to contribute code, patches or improvements to `COBREXA.jl`, please
read the basic [contribution guidelines and hints.](howToContribute.md).

## Acknowledgements

`COBREXA.jl` is developed at the Luxembourg Centre for Systems Biomedicine of
the University of Luxembourg ([uni.lu/lcsb](https://www.uni.lu/lcsb)),
cooperating with the Institute for Quantitative and Theoretical Biology at the Heinrich
Heine University in Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/)).

The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://www.permedcoe.eu/)) agreement no. 951773.

```@raw html
<img src="assets/cobrexa.svg" alt="COBREXA logo" style="height:8ex;width:auto">   <img src="assets/unilu.svg" alt="Uni.lu logo" style="height:8ex;width:auto">   <img src="assets/lcsb.svg" alt="LCSB logo" style="height:8ex;width:auto">   <img src="assets/hhu.svg" alt="HHU logo" style="height:8ex;width:auto">   <img src="assets/qtb.svg" alt="QTB logo" style="height:8ex;width:auto">   <img src="assets/permedcoe.svg" alt="PerMedCoE logo" style="height:8ex;width:auto">
```
