
# Modifying and saving the models

Making a small modification to the model and reanalyzing them is often a useful
way to explore how the constrains work together and what is the modification
freedom in the model.

With `COBREXA.jl`, you have 2 main choices of making such modifications:
- you can manually change the model structures
- many analysis functions accept special arguments that make the modifications
  in a declarative way and "on the fly", without requiring to manually interact
  with the model.

!!! tip "Notebook available"
    The available notebooks demonstrate
    [model export and serialization](../notebooks/1_loading_converting_saving.md)
    and various model modifications
    ([1](../notebooks/3_basic_stdmodel_usage.md),
    [2](../notebooks/4_basic_core_coupled_usage.md),
    [3](../notebooks/5_basic_stdmodel_construction.md)).

## Manual modifications

Certain model types, including [`CoreModel`](@ref) and [`StandardModel`](@ref),
are composed of mutable structures that you are free to modify as you want.

[`CoreModel`](@ref) consists of sparse matrices and vectors that describe the
model precisely. For example, modifying a bound of the reaction is as simple as
writing to the `.xl` or `.xu` (**l**ower and **u**pper bound for **x**) vectors
in the structure:
```
using COBREXA
m = load_model(CoreModel, "e_coli_core.xml")
m.xl[3] = 0.0
```

The available field names can be listed using e.g. `fieldnames(CoreModel)`, or
more conveniently by employing the completion in the shell:
```
julia> m.   # press [Tab]
S    b     c     mets  rxns  xl    xu
```

With `CoreModel`, you may need to find the proper metabolites by identifier.
For that, you may examine the [`reactions`](@ref) and [`metabolites`](@ref) of
the model, e.g. using
```
indexin(["M_nadh_c", "M_co2_e"], metabolites(m))
```
...which will return the numeric indexes of NADH and CO₂ metabolites. These can
be used e.g. to change the "balance" of the metabolites in the model:
```
m.b[64] = -1      # model will be losing 1 flux unit of CO₂
```
...or to modify existing reaction (here with index 5) directly in stoichiometry matrix:
```
m.S[5,8] = -1
m.S[5,64] = 1
```

While this works well if you are used to work with matrix-like representations
of the model, it is not really convenient if you want to change the reactions
and models in an easy way. [`StandardModel`](@ref) is structured in a much more
"object-oriented" way that makes the manual modifications easiler.

In particular, [`StandardModel`](@ref) consists of dictionaries of
[`Reaction`](@ref), [`Metabolite`](@ref) and [`Gene`](@ref) objects that may be
modified and indexed directly using their names. That way, the above
modifications may be done much more semantically, as follows:

```
m = load_model(StandardModel, "e_coli_core.xml")
m.reactions["R_TPI"].lb = 0.0                         # change lower bound of the reaction to 0
m.reactions["R_GLNS"].metabolites["M_nadh_c"] = -1.0  # update stoichiometry
m.reactions["R_GLNS"].metabolites["M_co2_e"] = 1.0
...
```

Finally, let's see the difference that the (biologically irrelevant) change has
done to the optimal flux of the model:

```
flux_balance_analysis_vec(m, GLPK.Optimizer) -
  flux_balance_analysis_vec(
    load_model(StandardModel, "e_coli_core.xml"), GLPK.Optimizer)
```

There are other functions that may be used to change the StandardModel in a
more systematic way. See the documentation of [`add!`](@ref), [`rm!`](@ref),
[`add_reaction!`](@ref), and [`set_bound`](@ref) for examples.

## Analysis modifiers

Some analysis functions, including [`flux_balance_analysis`](@ref) and
[`flux_variability_analysis`](@ref), accept a special argument `modifications`,
which is a list of descriptions of small changes that should be applied to the
model before modification.

These include e.g.:
- [`change_objective`](@ref) that sets a new optimization objective
- [`change_optimizer`](@ref) and [`change_optimizer_attribute`](@ref) that can
  be used to use a different `JuMP.jl` optimizer and set its parameters within
  the analysis
- [`change_constraint`](@ref) to change the flux bounds of a reaction
- [`knockout`](@ref) to disable reactions that depend on genes

This way, you can easily check out the model state when maximizing the rate of
"TALA" reaction:

```
m = load_model(StandardModel, "e_coli_core.xml")
flux_balance_analysis_dict(
  m, GLPK.Optimizer;
  modifications=[change_objective("R_TALA")])
```

...or try knocking out a gene combination that disables the "TALA" reaction
(see `m.reactions["R_TALA"].grr`):
```
flux_balance_analysis_dict(
  m, GLPK.Optimizer;
  modifications=[knockout(["G_b0008", "G_b2464"])])
```

...or do both at once; knock out some other genes and try maximizing the
reaction rate:
```
flux_balance_analysis_dict(
  m, GLPK.Optimizer;
  modifications=[
    knockout(["G_s0001"]),
    change_objective("R_TALA"),
  ])
```

## Exporting the modified models in native formats

Manually modified models can be exported in standard formats so that they can
be examined in other environments, or just made accessible for publication.

`COBREXA.jl` supports export of MATLAB-like and JSON models. Simply use
[`save_model`](@ref):

```
save_model(m, "myModel.json")
save_model(m, "myModel.mat")
```

The function automatically guesses the appropriate model type to write into the
file from extension; if required, you can choose the type of saved data
manually by using [`save_json_model`](@ref) and [`save_mat_model`](@ref).


## Using `Serialization` for quick&efficient model storage

If you save the model just "for yourself", for example just for the use in
later analysis, it may be inconvenient (and unnecessarily inefficient) to
encode and decode the models into&from the external format. Moreover, certain
model types (such as [`CoreModelCoupled`](@ref)) cannot be fully represented
in all model formats, thus increasing the chance for accidental data loss.

Instead of that, we recommend `Serialization` package that provides a
straightforward way to write out _any_ Julia data structure to the disk, using
a very efficient data format that can be written and read very quickly.

With any model in `m`, you can write it to disk as follows:
```
using Serialization
open(f -> serialize(f, m), "myModel", "w")
```
...and read it back with:
```
m = deserialize("myModel")
```

One great advantage of `Serialization` is speed -- models with millions of
reactions are usually loaded and saved with minimal overhead in less than a
second.

!!! note "Limits of `Serialization`"
    Serialized models are great for quickly exchanging data objects between
    analysis steps. The avoided need for re-encoding can save you a great deal
    of analysis time that can be used for better purposes. Despite of that, do
    not rely on the stability of the serialized format -- it often changes
    between Julia versions, and the data stored in one version may not open
    easily after an upgrade. In short, use serialized data within one workflow,
    and use standard and stable external formats for publishing and storing the
    data beyond the scope of a single analysis workflow.
