<br>
<div align="center">
    <img src="docs/src/assets/header.svg?maxAge=0" width="80%">
</div>

# COnstraint-Based Reconstruction and EXascale Analysis

[docs-img-stable]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url-stable]: https://lcsb-biocore.github.io/COBREXA.jl

[docs-img-dev]: https://img.shields.io/badge/docs-latest-0af.svg
[docs-url-dev]: https://lcsb-biocore.github.io/COBREXA.jl/dev/

[docker-url]: https://hub.docker.com/r/lcsbbiocore/cobrexa.jl
[docker-img]: https://img.shields.io/docker/image-size/lcsbbiocore/cobrexa.jl

[ci-img]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml/badge.svg?branch=master
[ci-url]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl/branch/master/graph/badge.svg?token=H3WSWOBD7L
[cov-url]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl

[contrib-img]: https://img.shields.io/badge/contributions-start%20here-green
[contrib-url]: https://lcsb-biocore.github.io/COBREXA.jl/stable/howToContribute/

[repostatus-url]: https://www.repostatus.org/#active
[repostatus-img]: https://www.repostatus.org/badges/latest/active.svg

| **Documentation** | **Tests** | **Coverage** | **How to contribute?** | **Project status** |
|:---:|:---:|:---:|:---:|:---:|
| [![docs-img-stable]][docs-url-stable] [![docs-img-dev]][docs-url-dev] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] | [![contrib][contrib-img]][contrib-url] | [![repostatus-img]][repostatus-url] |

This package provides constraint-based reconstruction and analysis tools for
exa-scale metabolic modeling in Julia.

## How to get started

### Prerequisites and requirements

- **Operating system**: Use Linux (Debian, Ubuntu or centOS), MacOS, or Windows
  10 as your operating system. `COBREXA` has been tested on these systems.
- **Julia language**: In order to use `COBREXA`, you need to install Julia 1.0
  or higher. Download and follow the installation instructions for Julia
  [here](https://julialang.org/downloads/).
- **Hardware requirements**: `COBREXA` runs on any hardware that can run Julia,
  and can easily use resources from multiple computers interconnected on a
  network. For processing large datasets, you are required to ensure that the
  total amount of available RAM on all involved computers is larger than the
  data size.
- **Optimization solvers**: `COBREXA` uses
  [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) to formulate optimization
  problems and is compatible with all [`JuMP` supported
  solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
  However, to perform analysis at least one of these solvers needs to be
  installed on your machine. For a pure Julia implementation, you may use e.g.
  [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl), but other solvers (GLPK,
  Gurobi, ...) work just as well.

:bulb: If you are new to Julia, it is advisable to [familiarize yourself with
the environment
first](https://docs.julialang.org/en/v1/manual/getting-started/).  Use the
Julia [documentation](https://docs.julialang.org) to solve various
language-related issues, and the [Julia package manager
docs](https://julialang.github.io/Pkg.jl/v1/getting-started/) to solve
installation-related difficulties. Of course, [the Julia
channel](https://discourse.julialang.org/) is another fast and easy way to find
answers to Julia specific questions.

### Quick start guide

<!--quickstart_begin-->
You can install COBREXA from Julia repositories. Start `julia`, **press `]`** to
switch to the Packaging environment, and type:
```
add COBREXA
```

You also need to install your favorite solver supported by `JuMP.jl`, typing
e.g.:
```
add Tulip
```

Alternatively, you may use [prebuilt Docker and Apptainer images](#prebuilt-images).

When the packages are installed, switch back to the "normal" julia shell by
pressing Backspace (the prompt should change color back to green). After that,
you can download [a SBML model from the
internet](http://bigg.ucsd.edu/models/e_coli_core) and perform a
flux balance analysis as follows:

```julia
using COBREXA   # loads the package
using Tulip     # loads the optimization solver

# download the model
download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

# open the SBML file and load the contents
model = load_model("e_coli_core.xml")

# run a FBA
fluxes = flux_balance_analysis_dict(model, Tulip.Optimizer)
```

The variable `fluxes` will now contain a dictionary of the computed optimal
flux of each reaction in the model:
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
  "R_G6PDH2r"     => 4.95999
  "R_TALA"        => 1.49698
  ⋮               => ⋮
```

#### Model variant processing

The main feature of COBREXA.jl is the ability to easily specify and process
many analyses in parallel. To demonstrate, let's see how the organism would perform if
some reactions were disabled independently:

```julia
# convert to a model type that is efficient to modify
m = convert(StandardModel, model)

# find the model objective value if oxygen or carbon dioxide transports are disabled
screen(m, # the base model
    variants=[ # this specifies how to generate the desired model variants
        [], # one with no modifications, i.e. the base case
        [with_changed_bound("R_O2t", lower=0.0, upper=0.0)], # disable oxygen
        [with_changed_bound("R_CO2t", lower=0.0, upper=0.0)], # disable CO2
        [with_changed_bound("R_O2t", lower=0.0, upper=0.0),
	        with_changed_bound("R_CO2t", lower=0.0, upper=0.0)], # disable both
    ],
    # this specifies what to do with the model variants (received as the argument `x`)
    analysis = x ->
        flux_balance_analysis_dict(x, Tulip.Optimizer)["R_BIOMASS_Ecoli_core_w_GAM"],
)
```
You should receive a result showing that missing oxygen transport makes the
biomass production much harder:
```julia
4-element Vector{Float64}:
 0.8739215022674809
 0.21166294973372796
 0.46166961413944896
 0.21114065173865457
```

Most importantly, such analyses can be easily specified by automatically
generating long lists of modifications to be applied to the model, and
parallelized.

Knocking out each reaction in the model is efficiently accomplished:

```julia
# load the task distribution package, add several worker nodes, and load
# COBREXA and the solver on the nodes
using Distributed
addprocs(4)
@everywhere using COBREXA, Tulip

# get a list of the workers
worker_list = workers()

# run the processing in parallel for many model variants
res = screen(m,
    variants=[
	# create one variant for each reaction in the model, with that reaction knocked out
        [with_changed_bound(reaction_id, lower=0.0, upper=0.0)]
	for reaction_id in reactions(m)
    ],
    analysis = model -> begin
	# we need to check if the optimizer even found a feasible solution,
	# which may not be the case if we knock out important reactions
    	sol = flux_balance_analysis_dict(model, Tulip.Optimizer)
	isnothing(sol) ? nothing : sol["R_BIOMASS_Ecoli_core_w_GAM"]
    end,
    # run the screening in parallel on all workers in the list
    workers = worker_list,
)
```

In result, you should get a long list of the biomass production for each
reaction knockout. Let's decorate it with reaction names:
```julia
Dict(reactions(m) .=> res)
```
...which should output an easily accessible dictionary with all the objective
values named, giving a quick overview of which reactions are critical for the
model organism to create biomass:
```julia
Dict{String, Union{Nothing, Float64}} with 95 entries:
  "R_ACALD"       => 0.873922
  "R_PTAr"        => 0.873922
  "R_ALCD2x"      => 0.873922
  "R_PDH"         => 0.796696
  "R_PYK"         => 0.864926
  "R_CO2t"        => 0.46167
  "R_EX_nh4_e"    => 1.44677e-15
  "R_MALt2_2"     => 0.873922
  "R_CS"          => 2.44779e-14
  "R_PGM"         => 1.04221e-15
  "R_TKT1"        => 0.864759
  ⋮             => ⋮
```
<!--quickstart_end-->

### Testing the installation

If you run a non-standard platform (e.g. a customized operating system), or if
you added any modifications to the `COBREXA` source code, you may want to run
the test suite to ensure that everything works as expected:

```julia
] test COBREXA
```

### Prebuilt images [![docker][docker-img]][docker-url]

Docker image is available from the docker hub as
[lcsbbiocore/cobrexa.jl][docker-url], and from GitHub container repository.
Download and use them as usual with docker:

```sh
docker run -ti --rm lcsbbiocore/cobrexa.jl:latest

# or alternatively from ghcr.io
docker run -ti --rm ghcr.io/lcsb-biocore/docker/cobrexa.jl:latest
```

In the container, you should get a `julia` shell with the important packages
already installed, and you may immediately continue the above tutorial from
`using COBREXA`.

Apptainer (aka Singularity) images are available from GitHub container
repository. To start one, run:
```sh
singularity run oras://ghcr.io/lcsb-biocore/apptainer/cobrexa.jl:latest
```
...which gives you a running Julia session with COBREXA.jl loaded.

If you require precise reproducibility, use a tag like `v1.2.2` instead of
`latest` (all releases since 1.2.2 are tagged this way).

<!--acknowledgements_begin-->
## Acknowledgements

`COBREXA.jl` is developed at the Luxembourg Centre for Systems Biomedicine of
the University of Luxembourg ([uni.lu/lcsb](https://wwwen.uni.lu/lcsb)),
cooperating with the Institute for Quantitative and Theoretical Biology at the Heinrich
Heine University in Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/)).

The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://permedcoe.eu/)) agreement no. 951773.
<!--acknowledgements_end-->

If you use COBREXA.jl and want to refer to it in your work, use the following citation format (also available as BibTeX in [cobrexa.bib](cobrexa.bib)):

> Miroslav Kratochvíl, Laurent Heirendt, St Elmo Wilken, Taneli Pusa, Sylvain Arreckx, Alberto Noronha, Marvin van Aalst, Venkata P Satagopam, Oliver Ebenhöh, Reinhard Schneider, Christophe Trefois, Wei Gu, *COBREXA.jl: constraint-based reconstruction and exascale analysis*, Bioinformatics, Volume 38, Issue 4, 15 February 2022, Pages 1171–1172, https://doi.org/10.1093/bioinformatics/btab782

<!--ack_logos_begin-->
<img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/cobrexa.svg" alt="COBREXA logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/unilu.svg" alt="Uni.lu logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/lcsb.svg" alt="LCSB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px" style="height:64px; width:auto">
<!--ack_logos_end-->
