# Model IO

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

rm(model_location) # hide
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.