# Model construction functions

## Functions for changing the models

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("reconstruction", file), readdir("../src/reconstruction"))
```

## Variant specifiers

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("reconstruction", "modifications", file), readdir("../src/reconstruction/modifications"))
```
