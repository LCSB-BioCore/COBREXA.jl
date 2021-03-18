# Model IO

## Reading constraint based models
Currently, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@example ioexample
using COBREXA

download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
model = read_model("e_coli_core.json")
rm("e_coli_core.json") # hide
model # pretty printing
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::CobraModel, file_location::String)
```

```@example ioexample
save_model(model, "e_coli_core2.json")
rm("e_coli_core2.json") # hide
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.
