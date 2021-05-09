
# Distributed processing and HPC environments

Distributed processing in Julia is represented mainly by the package
[`Distributed.jl`](https://docs.julialang.org/en/v1/stdlib/Distributed/).

`COBREXA.jl` is able to utilize this existing system to almost transparently
run the large parallelizable analyses on multiple CPU cores and multiple
computers connected through the network. Ultimately, the approach scales to
thousands of computing nodes in large HPC facilities.

Here, we give a short overview of how to work in the distributed environment
and utilize the resources for `COBREXA.jl` analyses.

## Starting the distributed workers

`COBREXA.jl` follows the structure imposed by the `Distributed` package: You
are operating a main (usually called "master") computation node, connect to
multiple other computers and start worker Julia processes there, and distribute
the workload across this community.

To start, you need to load the package and add a few processes. This starts 5
processes locally:

```
using Distributed
addprocs(5)
```

You may check that the workers are really there, using `workers()`. In this
case, it should give you a vector of _worker IDs_, very likely equal to
`[2,3,4,5,6]`.

If you have compute resources available via a network, you may connect these as
well, provided you have a secure shell (`ssh`) access to them. You will likely
want to establish a key-based authentication (refer to [ssh
documentation](https://www.openssh.com/manual.html)) to make the connection
easier.

With shell, check that you can `ssh` to a remote node and run Julia there:

```
user@pc> ssh server
...
user@server> julia
...
julia> _
```

If this works for you, you can add some workers that run on the `server` from
your Julia shell running on your `pc`. For example, this starts 20 workers on
the `server` and 10 workers on your friend's computer:

```
addprocs([('server', 20), ('joe_pc', 10)])
```

With this, you can schedule various computation on the workers; see the Julia
manual of [`Distributed`](https://docs.julialang.org/en/v1/stdlib/Distributed/)
for basic details. You may try various convenience packages, such as
[`DistributedArrays.jl`](https://github.com/JuliaParallel/DistributedArrays.jl)
and [`DistributedData.jl`](https://github.com/LCSB-BioCore/DistributedData.jl),
to process any data in a distributed fashion.

## Running a distributed analysis

While not all COBREXA functions may be parallelized naturally, these that do
will accept a special `workers` argument that specifies a list of worker IDs
where the computation should be distributed. For the value, you can use simply
`workers()`.

For example, [`flux_variability_analysis`](@ref) can naturally paralellize the
computation of all reactions's minima and maxima to finish the computation
faster. To enable the parallelization, you first need to make sure that all
workers have loaded both the COBREXA package and the optimizer:

```
using COBREXA, GLPK, Distributed
addprocs(10)                       # add any kind and number of processes here
@everywhere using COBREXA, GLPK    # loads the necessary packages on all workers
```

When the package is loaded and precompiled everywhere, you may load your model
and run the FVA with the `workers` parameter:

```
model = load_model("e_coli_core.xml")
result = flux_variability_analysis(model, GLPK.Optimizer; workers=workers())
```

With the extra computing capacity from `N` workers available, the FVA should be
computed roughly `N`-times faster.

!!! note "Distributed and parallel overhead"
    Communication of the workers with your Julia shell is not free. If the task
    that you are parallelizing is small and the model structure is very large,
    the distributed computation will actually spend most computation time just
    distributing the large model to the workers, and almost no time in
    executing the small parallel task. In such case, the performance will not
    improve by adding additional resources. You may want to check that the
    computation task is sufficiently large before investing the extra resources
    into the distributed execution.

## Interacting with HPC schedulers

Many researchers have access to institutional HPC facilities that allow
time-sharing of the capacity of a large computer cluster between many
researchers. Julia and `COBREXA.jl` work well within this environment; but your
programs require small additional customization to be able to find and utilize
the resources available from the HPC.

In our case, this reduces to a relatively complex task: You need to find out
how many resources were allocated for your task, and you need to add the remote
workers precisely at places that were allocated for your. Fortunately, the
package
[`ClusterManagers`](https://github.com/JuliaParallel/ClusterManagers.jl) can do
precisely that for us.

For simplicily, we will assume that your HPC is scheduled by
[Slurm](https://slurm.schedmd.com/).

Adding of the workers from Slurm is done as follows:
- you import the `ClusterManagers` package
- you find how many processes to spawn from the environment from `SLURM_NTASKS`
  environment variable
- you use the function `addprocs_slurm` to precisely connect to your allocated
  computational resources

The Julia script that does a parallel analysis may then start as follows:

```
using COBREXA, Distributed, ClusterManagers

available_workers = parse(Int, ENV["SLURM_NTASKS"])

addprocs_slurm(available_workers)

...
result = flux_variability_analysis(...; workers=workers())
...
```

After adding the Slurm workers, you may continue as if the workers were added
using normal `addprocs`, and (for example) run the
[`flux_variability_analysis`](@ref) as shown above.

!!! tip "What about the other HPC schedulers?"
    `ClusterManagers.jl` supports many other common HPC scheduling systems,
    including LFS, Sun Grid, SGE, PBS, and Scyld, in a way almost identical to
    Slurm. See the package documentation for details.

## Wrapping your script in a Slurm job

To be able to submit your script for later processing using the [`sbatch` Slurm
command](https://slurm.schedmd.com/sbatch.html), you need to wrap it in a small
"batch" script that tells Slurm how many resources the process needs.

Assuming you have a Julia computation script written down in `myJob.jl` and
saved on your HPC cluster's access node, the corresponding Slurm batch script
(let's call it `myJob.sbatch`) may look as follows:

```
#!/bin/bash -l
#SBATCH -n 100           # the job will require 100 individual workers
#SBATCH -c 1             # each worker will sit on a single CPU
#SBATCH -t 30            # the whole job will take less than 30 minutes
#SBATCH -J myJob         # the name of the job

module load lang/Julia   # this is usually required to make Julia available to your job

julia myJob.jl
```

To run the computation, simply run `sbatch myJob.sbatch` on the access node.
The job will be scheduled and eventually executed. You may watch `sacct` and
`squeue` in the meantime, to see the progress.

Remember that you need to explicitly save the result of your Julia script
computation to files, to be able to retrieve them later. Standard outputs of
the jobs are often mangled and discarded. If you still want to collect the
standard output, you may change the last line of the batch script to

```
julia myJob.jl > myJob.log
```

and collect the output from the log later. This is convenient especially if
logging various computation details using the `@info` and similar macros.
