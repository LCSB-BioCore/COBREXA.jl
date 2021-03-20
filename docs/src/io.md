# Model IO

## Reading constraint based models
Currently, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@example
using CobraTools

model_location = joinpath("..","..", "models", "iJO1366.json") 
model = read_model(model_location)
model # pretty printing
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::CobraTools.Model, file_location::String)
```

```@example
using CobraTools # hide
model_location = joinpath("..","..", "models", "iJO1366.json") # hide
model = read_model(model_location) # hide

# "e_coli_json_model.json" is the file name we are going to use to save the model
model_location = joinpath("e_coli_json_model.json")

# model is a CobraTools.Model object previously imported or created
save_model(model, model_location)

rm(model_location) # hide
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.