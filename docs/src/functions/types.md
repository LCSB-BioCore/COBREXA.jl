# Types

## Base types

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "types", "abstract", file), readdir("../src/base/types/abstract"))
```

## Model types and contents
```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "types", file), readdir("../src/base/types"))
```

## Model type wrappers
```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "types", "wrappers", file), readdir("../src/base/types/wrappers"))
```
