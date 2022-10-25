
# Writing custom optimizer modifications

Functions such as [`flux_balance_analysis`](@ref) internally create a JuMP
model out of the [`AbstractMetabolicModel`](@ref), and run the optimizer on that. To be
able to make some modifications on the JuMP model before the optimizer is
started, most of the functions accept a `modifications` argument, where one can
list callbacks that do the changes to the prepared optimization model.

The callbacks available in COBREXA.jl include functions that may help with
tuning the optimizer, or change the raw values in the linear model, such as:

- [`change_constraint`](@ref) and [`change_objective`](@ref)
- [`change_sense`](@ref), [`change_optimizer`](@ref), [`change_optimizer_attribute`](@ref)
- [`silence`](@ref)
- [`knockout`](@ref), [`add_crowding_constraints`](@ref)
- [`add_loopless_constraints`](@ref)

Compared to the [variant system](1_screen.md) and the [model
wrappers](4_wrappers.md), optimizer modifications are slightly more powerful
(they can do anything they want with the optimizer!), but do not compose well
-- it is very easy to break the semantics of the model or erase the previous
changes by carelessly adding the modifications.

Here, we show how to construct the modifications. Their semantics is similar to
the [variant-generating functions](1_screen.md), which receive a model (of type
[`AbstractMetabolicModel`](@ref)), and are expected to create another (modified) model.
Contrary to that, modifications receive both the [`AbstractMetabolicModel`](@ref) and a
JuMP model structure, and are expected to cause a side effect on the latter.

A trivial modification that does not do anything can thus be written as:

```julia
change_nothing() = (model, opt_model) -> println("Not touching anything.")
```

and applied as:
```julia
flux_balance_analysis(model, GLPK.Optimizer, modifications=[change_nothing()])
flux_variability_analysis(model, GLPK.Optimizer, modifications=[change_nothing()])
```

At the call time of the modifier function, `opt_model` is usually the model
that was returned from [`make_optimization_model`](@ref) -- refer to the
function for actual model layout. The function can freely change anything in
that model.

For demonstration, we show how to implement an impractical but illustrative
modification that adds an extra constraint that makes sure that all fluxes sum
to a certain value:

```julia
using JuMP

add_sum_constraint(total::Float64) =
    (model, opt_model) -> begin
        v = opt_model[:x]  # retrieve the variable vector
        @constraint(opt_model, total, sum(v) == total)  # create the constraint using JuMP macro
    end
```

The modification can be used at the expectable position:
```julia
v = flux_balance_analysis_vec(
    load_model("e_coli_core.xml"),
    GLPK.Optimizer,
    modifications = [add_sum_constraint(100.0)])

sum(v)   # should print ~100.0
```
