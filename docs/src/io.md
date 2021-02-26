# Model IO

## Model structure
Before reading or writing models, it is important to understand how they are represented internally.
Each model is a type of `CobraTools.Model`, which is composed of a model `id`, arrays of `Reaction`s, `Metabolite`s and `Gene`s, and a gene reaction rules (`grrs`) dictionary. 
The fields of `Reaction`, `Metabolite`, `Gene` types are shown below. 

When reading or writing, these fields are what is used by `CobraTools.jl`.
```@docs
Model
```
Note, the format of `grrs` in `CobraTools.Model` is a nested array, like [[g1, g2], [g3, g4], ...], indexed by a reaction `id` key. 
Each sub-array, e.g. [g1, g2], is composed of essential genes for the reaction to function. 
Thus, if rxn1 requires (g1 and g2) or (g3 and g4) to function, then this would be represented by rxn1 => [[g1, g2], [g3, g4]] in `grrs`. 
```@docs
Reaction
```
Note, the format of `grr` in `Reaction` is a string like `"(g1 and g2) or (g3 and g4)"`. 
This string is parsed into the format used by `grrs` in `CobraTools.Model`, as explained above.
Also note, the metabolites dictionary field of `Reaction` maps a `Metabolite` to its stoichiometrix coefficient.
```@docs
Metabolite
```
```@docs
Gene
```
FIY: `JuMP` also exports a `Model` type, so you need to qualify which `Model` you are referring to when making a new function.
## Reading constraint based models
Currently, SBML, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@setup modelread
using CobraTools

models_folder = joinpath(join(split(pwd(), "\\")[1:end-2], "\\"), "models")
model_location = joinpath(models_folder, "iJO1366.json") 
model = read_model(model_location) # import model so that the examples are simpler
```

```@example modelread
using CobraTools

# models_folder is a directory where models are stored
model_location = joinpath(models_folder, "iJO1366.json") 
model = read_model(model_location)
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::CobraTools.Model, file_location::String)
```

```@example modelread
using CobraTools

# "e_coli_json_model.json" is the file name we are going to use to save the model
model_location = joinpath(pwd(), "e_coli_json_model.json")

# model is a CobraTools.Model object previously imported or created
save_model(model, model_location)
```

## Problems
Please let me know when you run into model import/export problems by filing an issue on the repo.