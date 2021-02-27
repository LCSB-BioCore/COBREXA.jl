# CobraTools.jl
*CobraTools is a Julia package for constraint based reconstruction and analysis of metabolic models.*

[docs-img]:https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://stelmo.github.io/CobraTools.jl/dev
 
[ci-img]: https://github.com/stelmo/CobraTools.jl/actions/workflows/ci.yml/badge.svg?branch=master&event=push
[ci-url]: https://github.com/stelmo/CobraTools.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/stelmo/CobraTools.jl/branch/master/graph/badge.svg?token=3AE3ZDCJJG
[cov-url]: https://codecov.io/gh/stelmo/CobraTools.jl

[contrib]: https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat

[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE.md

[![][license-img]][license-url] [![contributions welcome][contrib]](https://github.com/stelmo/CobraTools.jl/issues)

| **Documentation** | **Tests** | **Coverage** |
|:--------------:|:-------:|:---------:|
| [![docs-img]][docs-url] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] |

This is package aims to provide constraint based reconstruction and analysis (COBRA) tools in the Julia environment, similar to Cobrapy in Python and the Cobra Toolbox in Matlab.
This package provides basic convenience functions, e.g. FBA, pFBA, sampling, model construction, etc.
More importantly, it also exposes the user to the core structures used in COBRA, e.g. the stoichiometric matrix, etc., so that custom optimization routines can be written as painlessly as possible (due in large part to JuMP). An alternative, [COBRA.jl](https://github.com/opencobra/COBRA.jl), but its scope is more restricted than `CobraTools.jl`.


## Installation

To install this package: `] add CobraTools`. See the documentation for more information.

## Quick Example

```julia
using CobraTools
using JuMP
using GLPK # pick any solver supported by JuMP

# import E. coli model
model = read_model("iJO1366.json") # models have pretty printing

# choose objective to maximize
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")

# Use convenience functions
sol = fba(model, biomass, GLPK.Optimizer; solver_attributes=Dict("msg_lev" => GLPK.GLP_MSG_OFF)) # classic flux balance analysis

# DIY
cbm, v, mb, ubs, lbs = build_cbm(model) # get the constraint based model (cbm) in JuMP format: S*v=b (mb: mass balance constraints) with lbs <= v <= ubs.
set_optimizer(cbm, GLPK.Optimizer) # use JuMP functions to set optimizer
set_optimizer_attribute(cbm, "msg_lev", GLPK.GLP_MSG_OFF) # use JuMP functions to set optimizer attributes
@objective(cbm, Max, v[model[biomass]]) # use index notation to get biomass equation index
optimize!(cbm)    
sol = map_fluxes(v, model) # map fluxes to reaction ids. 
```
More funcionality is described in the documention.

## Progress

- [x] Read JSON "Cobrapy" models
- [x] Read Matlab models
- [ ] Read SBML models
- [ ] Read YAML models
- [x] Write JSON models
- [x] Write Matlab models
- [ ] Write SBML models
- [ ] Write YAML
- [x] FBA
- [X] pFBA
- [ ] MOMA
- [ ] FVA
- [x] Implement sampling (hit and run)
- [x] Implement sampling (achr - kind of?)
- [ ] Single gene knockouts
- [ ] Double gene knockout
- [x] Equilibrator integration
- [x] Brenda integration (basic)
- [x] Reaction construction
- [x] Model modifications
- [ ] Distributed analysis (COBRA.jl integration?)

### Citations
1) Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 7(74). https://doi.org/10.1186/1752-0509-7-74
2) Heirendt, L., Arreckx, S., Pfau, T. et al. (2019). Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0. Nat Protoc 14, 639–702. https://doi.org/10.1038/s41596-018-0098-2
3) Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics, 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317
4) Chang, A., Jeske, L., Ulbrich, S., Hofmann, J., Koblitz, J., Schomburg, I., Neumann-Schaal, M., Jahn, D., Schomburg, D.. (2021). BRENDA, the ELIXIR core data resource in 2021: new developments and updates. Nucleic Acids Research, 49(D1). https://doi.org/10.1093/nar/gkaa1025
