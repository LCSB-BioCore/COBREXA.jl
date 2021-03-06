# CobraTools.jl
*CobraTools.jl is a Julia package for constraint based reconstruction and analysis of metabolic models.*

[docs-img]:https://img.shields.io/badge/docs-latest-blue.svg
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
This package provides basic convenience functions, e.g. model IO, construction, modification, and FBA, pFBA, sampling, etc. It can also be used to interface with [Equilibrator](http://equilibrator.weizmann.ac.il/) and [BRENDA](https://www.brenda-enzymes.org/).
More importantly, it also exposes the user to the core structures used in COBRA, e.g. the stoichiometric matrix, etc., so that custom optimization routines can be written as painlessly as possible (due in large part to [JuMP](https://jump.dev/)). An alternative Julia package hosted by openCOBRA, [COBRA.jl](https://github.com/opencobra/COBRA.jl), also exists, but its scope is more restricted than `CobraTools.jl`.


## Installation

To install this package: `] add https://github.com/stelmo/CobraTools.jl`. See the documentation for more information.

## Quick Example
Let's use `CobraTools.jl` to perform classic flux balance analysis on an *E. coli* genome-scale metabolic model.
```julia
using CobraTools
using JuMP
using Tulip # pick any solver supported by JuMP

# Import E. coli model (models have pretty printing)
model = read_model("iJO1366.json") 

# Choose objective to maximize (biomass is a reaction struct, which also has pretty printing)
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")

# FBA - use convenience function
sol = fba(model, biomass, Tulip.Optimizer)
```

If you are feeling more adventurous you can perform the optimization yourself using `JuMP`.
```julia
# Get the constraint based model (cbm) in JuMP format 
# cbm: S*v=b (mb: mass balance constraints) with lbs <= v <= ubs
cbm, v, mb, ubs, lbs = build_cbm(model)
# Use JuMP functions to optimize the constraint based model
set_optimizer(cbm, Tulip.Optimizer)
# Use index notation to get biomass equation index
@objective(cbm, Max, v[model[biomass]])
optimize!(cbm)    
# Map fluxes to reaction ids
sol = map_fluxes(v, model) 
```

If you are feeling even more adventurous you can do everything yourself!
```julia
# Get S*v = b with lbvec <= v <= ubvec from model
S, b, ubvec, lbvec = get_core_model(model) 
# Manually define the model in JuMP
cbm = JuMP.Model(Tulip.Optimizer)
nvars = size(S, 2)
v = @variable(cbmodel, v[1:nvars]) 
mb = @constraint(cbmodel, mb, S*v .== b) 
lbs = @constraint(cbmodel, lbs, lbs .<= v) 
ubs = @constraint(cbmodel, ubs, v .<= ubs) 
@objective(cbm, Max, v[model[biomass]])
optimize!(cbm)
sol = map_fluxes(v, model)
```
More funcionality is described in the documention, e.g. model construction and analysis in pure Julia.

## Progress

- [x] Read JSON "Cobrapy" models
- [x] Read Matlab models
- [ ] Read SBML models
- [x] Write JSON models
- [x] Write Matlab models
- [ ] Write SBML models
- [x] FBA
- [X] pFBA
- [ ] FVA
- [x] Implement sampling (hit and run)
- [x] Implement sampling (achr - kind of, fix constraint issue...)
- [ ] Single gene knockouts
- [ ] Double gene knockout
- [x] Equilibrator integration
- [x] Brenda integration (basic)
- [x] Reaction construction
- [x] Model construction
- [x] Model modifications
- [ ] Thermodynamic FBA (and related functions)



### Citations
1) Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 7(74). https://doi.org/10.1186/1752-0509-7-74
2) Heirendt, L., Arreckx, S., Pfau, T. et al. (2019). Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0. Nat Protoc 14, 639–702. https://doi.org/10.1038/s41596-018-0098-2
3) Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics, 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317
4) Chang, A., Jeske, L., Ulbrich, S., Hofmann, J., Koblitz, J., Schomburg, I., Neumann-Schaal, M., Jahn, D., Schomburg, D.. (2021). BRENDA, the ELIXIR core data resource in 2021: new developments and updates. Nucleic Acids Research, 49(D1). https://doi.org/10.1093/nar/gkaa1025
