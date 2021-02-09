# CobraTools.jl

This is package aims to provide basic constraint based reconstruction analysis (COBRA) tools in the Julia environment.

## Installation

To import Matlab formatted models you will need `MATLAB.jl` installed and a working Matlab. Refer to the documentation of MATLAB.jl for installation instructions.

The optimization solvers are implemented through `JuMP` and thus this package should be solver agnostic. All tests are conducted using `Gurobi.jl` but other solvers should work. 

To install this package: `] add CobraTools`.

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
- [ ] Implement sampling
- [ ] Single gene knockouts
- [ ] Double gene knockout
- [ ] Model construction tools
- [ ] Gibbs integration
- [ ] Distributed analysis (COBRA.jl integration?)

## Usage

### Read and writing models

### Basic analysis
FBA

FVA

pFBA

Knock-outs

### Sampling

### Gibbs 
The standard (temperature 25 °C, pressue 1 bar, concentration 1M) Gibbs energies of biochemical reactions at various pH levels (ionic strength = 0.1M) mapped to the KEGG database are made available from [Equilibrator](http://equilibrator.weizmann.ac.il/download). The raw .csv files downloaded from Equilibrator are located in `data`. Temperature and reactant concentration adjustments are made by ΔᵣG = ΔᵣG⁰ + RTln(Q).   

### Citations
Cobrapy

    Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 7(74). https://doi.org/10.1186/1752-0509-7-74

ACHR

CHRR

Equilibrator

    Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics, 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317

MetaNetX
    Bernard, T., Bridge A., Morgat A., Moretti, S., Xenarios, I., Pagni, M., (2014) Reconciliation of metabolites and biochemical reactions for metabolic networks, Briefings in Bioinformatics, 15(1), 123–135, https://doi.org/10.1093/bib/bbs058