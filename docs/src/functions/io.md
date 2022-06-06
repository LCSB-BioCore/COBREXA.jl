# Input and output

## File I/O and serialization

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("io", file), readdir("../src/io"))
```

## Pretty printing

```@autodocs
Modules = [COBREXA]
Pages = map(file -> joinpath("io", "show", file), readdir("../src/io/show"))
```
