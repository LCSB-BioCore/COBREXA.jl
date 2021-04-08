# Model IO

## Reading constraint based models
Currently, JSON and Matlab formatted models can be imported.

```@docs
read_model(file_location::String)
```

```@example ioexample
using COBREXA

if !isfile("e_coli_core.json")
  download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
end

model = read_model("e_coli_core.json")
model # pretty print the model
```

## Writing constraint based models
Currently, JSON and Matlab models can be exported.

```@docs
save_model(model::StandardModel, file_location::String)
```

```@example ioexample
save_model(model, "e_coli_core2.json")
```

## IO Problems?
Please let me know when you run into model import/export problems by filing an issue.
