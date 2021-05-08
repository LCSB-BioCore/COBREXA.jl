
# Loading and converting model data

COBREXA.jl supports several constraint-based model formats that can be loaded
with built-in functions. You can load SBML models that are compatible with
[libsbml](http://sbml.org/Software/libSBML), JSON models (such as the ones from
[CobraPy](https://github.com/opencobra/cobrapy)), and MATLAB-style models (such
as those from [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)).

These formats are commonly available from many model repositories, such as from
BIGG, as seen e.g. on the BIGG entry for the core [*E. Coli*
entry](http://bigg.ucsd.edu/models/e_coli_core). Here, we show how to load the
basic formats and work with such models.

!!! tip "Notebook available!"
    Example code for this tutorial is [available
    here](../notebooks/1_loading_converting_saving.md).

## Loading models from files

For most purposes, you should be able to open and load any model with
[`load_model`](@ref), which detects the file type from the extension (`.xml`,
`.json` and `.mat`), and calls the appropriate loading function. After loading
the `COBREXA.jl` library with `using COBREXA` and you may load the downloaded
model data into Julia as follows:

```
my_model = load_model("e_coli_core.xml")
```

You should see some information about the loaded model, possibly looking like
this:
```
Metabolic model of type JSONModel

  [9 ,  1]  =  1.0
  [51,  1]  =  1.0
  ⋮
  [57, 95]  =  1.0
  [59, 95]  =  -1.0
Number of reactions: 95
Number of metabolites: 72
```

If the file type can not be guessed from the file extension, use any of the
specific loader functions:

- [`load_sbml_model`](@ref) for SBML
- [`load_json_model`](@ref) for JSON
- [`load_mat_model`](@ref) for MATLAB models

All formats may store slightly different information. By default, COBREXA
attempts not to discard any information unless a conversion to a more strict
format is required. For example, the [`JSONModel`](@ref) (which is returned by
[`load_json_model`](@ref)) still holds the original JSON structure that you can
freely access for any specific purposes:

```
jm = load_json_model("e_coli_core.json")
jm.json["reactions"][1]
```

That should print out the first reaction in the model in a JSON-style scheme,
in our case the process catalyzed by phosphofructokinase:

```
Dict{String,Any} with 9 entries:
  "name"               => "Phosphofructokinase"
  "metabolites"        => Dict{String,Any}("adp_c"=>1.0,"atp_c"=>-1.0,"f6p_c"=>…
  "lower_bound"        => 0.0
  "id"                 => "PFK"
  "notes"              => Dict{String,Any}("original_bigg_ids"=>Any["PFK"])
  "gene_reaction_rule" => "b3916 or b1723"
  "upper_bound"        => 1000.0
  "subsystem"          => "Glycolysis/Gluconeogenesis"
  "annotation"         => Dict{String,Any}("ec-code"=>Any["2.7.1.11"],"metanetx…
```

[`MATModel`](@ref) and [`SBMLModel`](@ref) (returned by the respective loading
functions) contain similar "full" model information -- you can access the whole
MATLAB and SBML data and build on them without any restrictions.

## Converting to other model types

Despite JSON and SBML are great for storing and exchanging the models, the data
representation is not very suitable for analyzing the model and processing it
mathematically.

COBREXA.jl contains several model types that are much better suited  for
supporting the analysis tasks. You can use the following:

- [`CoreModel`](@ref), which represents the "core" of the optimization problem
  and the corresponding linear programming problem -- a sparse representation
  of the stoichiometric matrix, flux bounds vectors, objective vector, etc.
- [`StandardModel`](@ref) (a "standard" for COBREXA.jl), which represents a
  highly flexible, object-like, dictionary-based representation of a model that
  contains individual [`Reaction`](@ref)s, [`Metabolite`](@ref)s,
  [`Gene`](@ref)s, and other things.

!!! note "Conversion limitations and possible data loss"
    Because of the specifics of the format of each model structure, the
    conversion is not always able to preserve all information from the source
    data. You may need to check if any complicated and less-standard
    annotations are still present. If you require them, and either use a more
    complicated model, or collect them manually.

A loaded model can be converted to any other model type using the standard
Julia conversion:

```
cm = convert(CoreModel, jm)
```

You can also use a shortcut in [`load_model`](@ref) to convert the model to the
desired format in one command:

```
cm = load_model(CoreModel, "e_coli_core.xml")
```

With [`CoreModel`](@ref), the information is easily accessible in matrix form.
For example, `cm.S` now contains the sparse stoichiometric matrix, which you
can convert to a dense matrix and manipulate it in Julia as any other matrix:

```
Matrix(cm.S)
```
...should show you the (relatively empty) stoichiometry of the model.

[`StandardModel`](@ref) is more suitable for fine-grained access to individual
items of the model, perhaps closer to the SBML-style models. For example, you
can view and set reaction bounds as follows:

```
sm = load_model(StandardModel, "e_coli_core.json")
sm.reactions["PGI"].ub
```
...this prints the upper bound of the reaction (in this case, `1000.0`); you
can change it the usual way:
```
sm.reactions["PGI"].ub = 500
```

This change will naturally project to future analysis results.
