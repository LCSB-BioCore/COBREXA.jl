# Utilities

## Helper functions

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "utils", file), readdir("../src/base/utils"))
```

## Macro-generated functions and internal helpers

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "macros", file), readdir("../src/base/macros"))
```

## Logging and debugging helpers

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "logging", file), readdir("../src/base/logging"))
```
