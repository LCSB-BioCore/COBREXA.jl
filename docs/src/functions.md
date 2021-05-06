# Functions

## Base Types

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "types", "abstract", file), readdir("../src/base/types/abstract"))
```

## Model types and contents

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "types", file), readdir("../src/base/types"))
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

### Pretty printing

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("io", "show", file), readdir("../src/io/show"))
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

### Analysis modifications

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("reconstruction", "modifications", file), readdir("../src/reconstruction/modifications"))
```

### Flux sampling

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("analysis", "sampling", file), readdir("../src/analysis/sampling"))
```

## Miscellaneous utilities

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "utils", file), readdir("../src/base/utils"))
```

### Logging and debugging helpers

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("base", "logging", file), readdir("../src/base/logging"))
```
