# Functions

## Types

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("types", file), readdir("../src/types"))
```

## Base functions

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", file), readdir("../src/base"))
```

## File I/O and serialization

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("io", file), readdir("../src/io"))
```

## Model reconstruction

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("reconstruction", file), readdir("../src/reconstruction"))
```

## Analysis functions

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("analysis", file), readdir("../src/analysis"))
```
