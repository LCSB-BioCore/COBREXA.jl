# Model IO

## Reading constraint based models
Currently, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@example io
using COBREXA

# Download a model from the BIGG database
download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Read that model into Julia
model = read_model("e_coli_core.json")
model # pretty printing
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::CobraModel, file_location::String)
```

```@example io
save_model(model, "e_coli_core_saved.json")
rm("e_coli_core2.json") # hide
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.
