<br>
<div align="center">
    <img src="docs/src/assets/header.svg?maxAge=0" width="80%">
</div>

# COnstraint-Based Reconstruction and EXascale Analysis

[docs-img-stable]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url-stable]: https://lcsb-biocore.github.io/COBREXA.jl

[docs-img-dev]: https://img.shields.io/badge/docs-latest-0af.svg
[docs-url-dev]: https://lcsb-biocore.github.io/COBREXA.jl/dev/

[docs-url-quickstart]: https://lcsb-biocore.github.io/COBREXA.jl/stable/quickstart/
[docs-url-examples]: https://lcsb-biocore.github.io/COBREXA.jl/stable/examples/

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

**This repository contains the "legacy" version of COBREXA (version 1.x).
COBREXA development is continuing at
[github.com/COBREXA](https://github.com/COBREXA/) with a much-improved version
2.x, and many other new packages. Unless you need version 1.x for compatibility
reasons, we recommend switching to version 2!**

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

### Quick start

[COBREXA.jl documentation][docs-url-stable]
is available online (also for
[development version][docs-url-dev]
of the package).

<!--quickstart_begin-->
You can install COBREXA from Julia repositories. Start `julia`, **press `]`** to
switch to the Packaging environment, and type:
```
add COBREXA@1
```

(Note: the above command installs the "1.x" version of COBREXA, which is what
is described in this repository and the associated documentation. If you are
not limited by compatibility and similar reasons, consider updating to [COBREXA
2](https://github.com/COBREXA/COBREXA.jl).)

You also need to install your favorite solver supported by `JuMP.jl` (such as
Gurobi, Mosek, CPLEX, GLPK, Clarabel, etc., see a [list
here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)).  For
example, you can install `Tulip.jl` solver by typing:
```
add Tulip
```

Alternatively, you may use [prebuilt Docker and Apptainer images](#prebuilt-images).

If you are running COBREXA.jl for the first time, it is very likely that upon
installing and importing the packages, your Julia installation will need to
precompile their source code from the scratch. In fresh installations, the
precompilation process should take less than 5 minutes.

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
<!--quickstart_end-->

The main feature of COBREXA.jl is the ability to easily specify and process a
huge number of analyses in parallel. You. You can have a look at a
[longer guide that describes the parallelization and screening functionality][docs-url-quickstart],
or dive into the [example analysis workflows][docs-url-examples].

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

If you use COBREXA.jl and want to refer to it in your work, use the following
citation format (also available as BibTeX in [cobrexa.bib](cobrexa.bib)):

> Miroslav Kratochvíl, Laurent Heirendt, St Elmo Wilken, Taneli Pusa, Sylvain Arreckx, Alberto Noronha, Marvin van Aalst, Venkata P Satagopam, Oliver Ebenhöh, Reinhard Schneider, Christophe Trefois, Wei Gu, *COBREXA.jl: constraint-based reconstruction and exascale analysis*, Bioinformatics, Volume 38, Issue 4, 15 February 2022, Pages 1171–1172, https://doi.org/10.1093/bioinformatics/btab782

<!--ack_logos_begin-->
<img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/cobrexa.svg" alt="COBREXA logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/unilu.svg" alt="Uni.lu logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/lcsb.svg" alt="LCSB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/COBREXA.jl/stable/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px" style="height:64px; width:auto">
<!--ack_logos_end-->
