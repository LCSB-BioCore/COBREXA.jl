# Model Structure
Before reading or writing models, it is important to understand how they are represented internally.
Each model is a type of `CobraTools.Model`, which is composed of a model `id`, and arrays of `Reaction`s, `Metabolite`s and `Gene`s. 
The fields of `Reaction`, `Metabolite`, `Gene` types are shown below. 

When reading or writing, these fields are what is used by `CobraTools.jl`.
```@docs
Model
```
```@docs
Reaction
```
Note, the format of `grr` in `Reaction` should be a nested array, like [[g1, g2], [g3, g4], ...]. 
Each sub-array, e.g. [g1, g2], is composed of essential genes (`g1::CobraTools.Gene`, etc.) for the reaction to function. 
Thus, if the reaction requires (g1 and g2) or (g3 and g4) to function, then this would be represented by [[g1, g2], [g3, g4]] in `grr`. 

Also note, the metabolites dictionary field of `Reaction` maps a `Metabolite` to its stoichiometrix coefficient.
```@docs
Metabolite
```
```@docs
Gene
```
FIY: `JuMP` also exports a `Model` type, so you need to qualify which `Model` you are referring to when making a new function.