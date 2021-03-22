# Model Structure

Before reading, writing, or building models, it is important to understand how
they are represented internally in `COBREXA.jl`.

Each model is a struct of the type `CobraModel`, which is composed of a model
`id`, and arrays of `Reaction`s, `Metabolite`s and `Gene`s.

```@docs
CobraModel
```

The fields of `Reaction`, `Metabolite`, `Gene` types are shown below.  When
reading, writing, building or analysing models, these fields are what is used
by `COBREXA.jl`.

```@docs
Reaction
```

Note, the format of `grr` (gene reaction rule) in `Reaction` should be a nested
array, like `[[g1, g2], [g3, g4], ...]`.

Each sub-array, e.g. `[g1, g2]`, is composed of essential genes
(`g1::COBREXA.Gene`, etc.) for the reaction to function.

Thus, if the reaction requires (`g1` and `g2`) or (`g3` and `g4`) to function,
then this would be represented by `[[g1, g2], [g3, g4]]` in `grr`. Also note,
the metabolites dictionary field of `Reaction` maps a `Metabolite` to its
stoichiometrix coefficient, i.e. it represents the reaction equation.

Also note, the format used the `annotation` field in `Reaction`, `Metabolite`
and `Gene` should always be a dictionary with string ids mapped to vectors of
entries EXCEPT for "sbo" terms, which are not vectors but strings, e.g.:
```
reaction.annotation = Dict(
  "bigg.reaction" => ["PFK"],
  "rhea" => ["16111", "16109", "16110", "16112"],
  "sbo" => "SBO:0000176",
  "ec-code" => ["2.7.1.11"],
)
```
...is correct, but any other format will cause issues.

```@docs
Metabolite
```

```@docs
Gene
```
