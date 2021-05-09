```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# Constraint-Based Reconstruction and EXascale Analysis

COBREXA is a toolkit for working with large constraint-based metabolic models,
and running very large number of analysis tasks on these models in parallel.
Its main purpose is to make the methods of Constraint-based Reconstruction and
Analysis (COBRA) scale to problem sizes that require the use of huge computer
clusters and HPC environment, which allows to realistically approach the
pre-exascale-sized models.

In the package, you will find the usual COBRA-like functions that inteface to
the underlying linear programming solvers. We use
[`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface for
many solvers; you can plug in whatever compatible solver you want, including
the popular `Tulip.jl`, `GLPK.jl`, `OSQP.jl`, and `Gurobi`.

## Quick start guide

<!--insert_quickstart-->

## Tutorials

Detailed tutorial contents is [available here](tutorials.md).

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

If you aim to contribute code, patches or improvements to `COBREXA.jl` read the
basic guidelines at a separate page with [contribution guidelines and
hints](howToContribute.md).
