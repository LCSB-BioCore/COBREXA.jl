<div align="center">
    <img src="docs/src/assets/header.svg?maxAge=0" width="80%">
</div>
<br>

# Constraint-Based Reconstruction and EXascale Analysis

[docs-img]:https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: http://lcsb-biocore.github.io/COBREXA.jl

[ci-img]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml/badge.svg?branch=master&event=push
[ci-url]: https://github.com/LCSB-BioCore/COBREXA.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl/branch/master/graph/badge.svg?token=3AE3ZDCJJG
[cov-url]: https://codecov.io/gh/LCSB-BioCore/COBREXA.jl

[contrib]: https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat

[![contributions welcome][contrib]](https://github.com/LCSB-BioCore/COBREXA.jl/issues)

| **Documentation** | **Tests** | **Coverage** |
|:--------------:|:-------:|:---------:|
| [![docs-img]][docs-url] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] |

This is package aims to provide constraint based reconstruction and analysis tools at the exa-scale in Julia.


## Installation

To install this package: `] add ???`. See the documentation for more information.

## Quick Example
Let's use `COBREXA.jl` to perform classic flux balance analysis on an *E. coli* community.
```julia
using COBREXA
using JuMP
using Tulip # pick any solver supported by JuMP

# Import E. coli models (models have pretty printing)
model_1 = read_model("iJO1366.json")
model_2 = read_model("iJO1366.json")
model_3 = read_model("iJO1366.json")

# Build an exascale model
exascale_model = join(model_1, model_2, model_3,...)
```
More funcionality is described in the documention, e.g. model construction and analysis in pure Julia.

### Citations
1) Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 7(74). https://doi.org/10.1186/1752-0509-7-74
2) Heirendt, L., Arreckx, S., Pfau, T. et al. (2019). Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0. Nat Protoc 14, 639–702. https://doi.org/10.1038/s41596-018-0098-2
3) Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics, 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317
4) Chang, A., Jeske, L., Ulbrich, S., Hofmann, J., Koblitz, J., Schomburg, I., Neumann-Schaal, M., Jahn, D., Schomburg, D.. (2021). BRENDA, the ELIXIR core data resource in 2021: new developments and updates. Nucleic Acids Research, 49(D1). https://doi.org/10.1093/nar/gkaa1025
