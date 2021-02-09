# CobraTools.jl

This is package aims to provide basic constraint based reconstruction analysis (COBRA) tools in the Julia environment.

## Installation
To import SBML formatted models you will need both `PyCall` and `python-libsbml` installed. See [here](https://stochasticreveller.wordpress.com/2016/08/02/sbml-and-the-julia-programming-language/) for some installation hints.

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

    Ebrahim, A., Lerman, J.A., Palsson, B.O. & Hyduke, D. R. (2013) COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Syst Biol 7, 74. https://doi.org/10.1186/1752-0509-7-74

ACHR

CHRR

Equilibrator

    Noor, E., Bar-Even, A., Flamholz, A., Lubling, Y., Davidi, D., & Milo, R. (2012). An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics (Oxford, England), 28(15), 2037–2044. https://doi.org/10.1093/bioinformatics/bts317