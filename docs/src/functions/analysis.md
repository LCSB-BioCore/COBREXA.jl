# Analysis functions

## Common analysis functions

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("analysis", file), readdir("../src/analysis"))
```

## Sampling

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("analysis", "sampling", file), readdir("../src/analysis/sampling"))
```

## Analysis modifiers

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("analysis", "modifications", file), readdir("../src/analysis/modifications"))
```
