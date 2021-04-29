# Tutorials

## Overview
TODO

## Tutorial 1: Loading models
COBREXA can load models stored in `.mat`, `.json`, and `.xml` formats (with the latter
denoting SBML formatted models). Each model type is read into memory as a `MATModel`, 
`JSONModel` and `SBMLModel` respectively, with the model type being inferred from the file extension.
```@example
mat_model = load_model("e_coli_core.mat")
json_model = load_model("e_coli_core.json")
sbml_model = load_model("e_coli_core.xml")
```
Often it is convenient to interconvert models between different types, e.g. `MATModel` -> `SBMLModel`.
COBREXA makes use of a generic interface to faciliate this. This makes it possible to interconvert
between all model types, however this conversion is not guaranteed to be lossless. In all cases
the most

## Tutorial 2: Joining and reconstructing models
TODO

## Tutorial 3: Flux balance analysis
TODO

## Tutorial 4: Flux variability analysis
TODO

## Tutorial 5: Parsimonious flux balance analysis
TODO

## Tutorial 6: Loopless flux balance analysis
TODO

## Tutorial 7: Sampling
TODO
