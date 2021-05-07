```@raw html
<br>
<div align="center">
    <img src="assets/header.svg?maxAge=0" width="80%">
</div>
<br>
```

# COnstraint-Based Reconstruction and EXascale Analysis

COBREXA is a toolkit for working with large constraint-based metabolic models,
and running a very large number of analysis tasks on these models in parallel.
Its main purpose is to make the methods of Constraint-based Reconstruction and
Analysis (COBRA) scale to problem sizes that require the use of huge computer
clusters and HPC environments. This allows these analysis techniques to
realistically approach the pre-exascale-sized models.

In the package, you will find the usual COBRA-like functions that interface with
the underlying optimizers. We use
[`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) as the unified interface to
manipulate these optimizers. Thus, you can plug in any compatible optimizer
you want, including the popular `Tulip.jl`, `GLPK.jl`, `OSQP.jl`, and `Gurobi`
optimizers, and everything will just work.

## Quick start guide

You can install COBREXA using the Julia package manager. Start `julia`, **press `]`** to
switch to the Packaging environment, and type:
```
add COBREXA
```

You may also need to install your favorite solver, typing e.g.:
```
add GLPK
```

When the packages are installed, switch back to the "normal" julia shell by
pressing backspace (the prompt should change color back to green). After that,
you can download [an SBML model from the
internet](http://bigg.ucsd.edu/models/e_coli_core) and perform
flux balance analysis as follows:

```julia
using COBREXA   # loads the package
using GLPK      # loads your favorite solver

isfile("e_coli_core.xml") || download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

model = load_model("e_coli_core.xml")

fluxes = flux_balance_analysis_dict(model, GLPK.Optimizer)
```

The variable `fluxes` will now contain a dictionary of the computed
flux distribution maximizing the biomass objective function of each reaction in the model:
```
Dict{String,Float64} with 95 entries:
  "R_EX_fum_e"    => 0.0
  "R_ACONTb"      => 6.00725
  "R_TPI"         => 7.47738
  "R_SUCOAS"      => -5.06438
  "R_GLNS"        => 0.223462
  "R_EX_pi_e"     => -3.2149
  "R_PPC"         => 2.50431
  "R_O2t"         => 21.7995
  "R_G6PDH2r"     => 4.95998
  "R_TALA"        => 1.49698
  ⋮               => ⋮
```

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

If you want to contribute, please read [the contribution guidelines and hints](howToContribute.md).
