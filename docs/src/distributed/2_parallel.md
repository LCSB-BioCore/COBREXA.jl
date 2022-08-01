
# Local parallel processing

To run an analysis in parallel, you first need to load the `Distributed` package and add a few worker processes. For example, you may start 5 local processes (that may utilize 5 CPUs) like this:

```julia
using Distributed
addprocs(5)
```

!!! note "`Distributed.jl` installation"
    `Distributed.jl` usually comes pre-installed with Julia distribution, but
    you may still need to "enable" it by typing `] add Distributed`.

You may check that the workers are really there, using `workers()`. In this
case, it should give you a vector of _worker IDs_, very likely equal to
`[2,3,4,5,6]`.

Each of the processes contains a self-sufficient image of Julia that can act
independently; in turn the additional processes also consume some memory. Each
process with loaded `COBREXA.jl` and a simple solver such as GLPK may consume
around 500MB of RAM, which should be taken into account when planning the
analysis scale.

Packages (COBREXA and your selected solver) must be loaded at all processes,
which you can ensure using the "everywhere" macro:
```julia
@everywhere using COBREXA, GLPK
```

Utilizing the prepared worker processes is then straightforward: You pass the
list of workers to the selected analysis function using the `workers` keyword
argument, and the parallel processing is automatically orchestrated for you.

```julia
model = load_model("e_coli_core.xml")
result = flux_variability_analysis(model, GLPK.Optimizer; workers=workers())
```
