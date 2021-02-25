# CobraTools.jl
*CobraTools is a Julia package for constraint based reconstruction and analysis of metabolic models.*

This is package aims to provide constraint based reconstruction and analysis (COBRA) tools in the Julia environment.
This package provides basic convenience functions, e.g. FBA, pFBA, sampling, model construction, etc.
More importantly, it also exposes the user to the core machinery used in this type of analysis, e.g. the stoichiometric matrix, so that custom optimization routines can be written as painlessly as possible (due in large part to JuMP). 

## Installation

To install this package: `] add CobraTools`. See the documentation for more information.

## Quick Example

TODO.

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
- [x] Gibbs integration
- [x] Brenda integration (basic)
- [x] Reaction construction
- [ ] Model modifications
- [ ] Distributed analysis (COBRA.jl integration?)

### Citations
1) Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 7(74). https://doi.org/10.1186/1752-0509-7-74
2) Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics, 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317