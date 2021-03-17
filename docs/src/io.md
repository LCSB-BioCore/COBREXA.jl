# Model IO

## Reading constraint based models
Currently, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@example ioexample
using COBREXA

model_location = download("http://bigg.ucsd.edu/static/models/iJO1366.json" ,"iJO1366.json") 
model = read_model(model_location)
model # pretty printing
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::CobraModel, file_location::String)
```

```@example ioexample
model_location = joinpath("e_coli_json_model.json")
save_model(model, model_location)
rm(model_location) # hide
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.