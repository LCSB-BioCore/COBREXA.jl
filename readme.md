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
- [ ] FBA
- [ ] pFBA
- [ ] MOMA
- [ ] FVA
- [ ] Implement sampling
- [ ] Single gene knockouts
- [ ] Double gene knockout
- [ ] Model construction tools
- [ ] Distributed analysis (COBRA.jl integration?)

## Usage

### Read and writing models

### Basic analysis
FBA

FVA

pFBA

Knock-outs

### Sampling