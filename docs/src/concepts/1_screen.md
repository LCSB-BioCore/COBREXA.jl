
# Screening many model variants

A major goal of COBREXA.jl is to make exploring of many model variants easy and
fast.

One main concept that can be utilized for doing that is implemented in the
function [`screen`](@ref), which takes your model, a list of model _variants_
that you want to explore by some specified _analysis_, and schedules the
analysis of the model variants parallelly on the available distributed workers.

In its most basic form, the "screening" may use the slightly simplified variant
of [`screen`](@ref) that is called [`screen_variants`](@ref), which works as
follows:

```julia
m = load_model(ObjectModel, "e_coli_core.json")

screen_variants(
    m,    # the model for screening
    [
        [],    # a variant with no modifications
        [with_changed_bound("CO2t", lb = 0, ub = 0)],  # disable CO2 transport
        [with_changed_bound("O2t", lb = 0, ub = 0)],  # disable O2 transport
        [with_changed_bound("CO2t", lb = 0, ub = 0), with_changed_bound("O2t", lb = 0, ub = 0)],  # disable both transports
    ],
    m -> flux_balance_analysis_dict(m, Tulip.Optimizer)["BIOMASS_Ecoli_core_w_GAM"],
)
```
The call specifies a model (the `m` that we have loaded) that is being tested,
then a vector of model variants to be created and tested, and then the analysis
that is being run on each variant -- in this case, we find an optimal steady
state of each of the variants, and check out the biomass production rate at
that state. In this particular case, we are checking what will be the effect of
disabling combinations of CO2 transport and O2 transport in the cells. For
that, we get the following result:
```
4-element Vector{Float64}:
 0.8739215022678488
 0.46166961413944896
 0.21166294973372135
 0.21114065173865518
```

The numbers are the biomass production rates for the specified variants. We can
see that disabling O2 transport really does not help the organism much.

## Variant specification

In the above example, we have specified 4 variants, thus the analysis returned
4 different results that correspond with the specifications. Let us have a look
at the precise format of the specification and result.

Importantly, the `variants` argument is of type `Array{Vector{Any}}`, meaning
that it can be an array of any dimensionality that contains vectors. Each of the
vectors specifies precisely one variant, possibly with more modifications
applied to the model in sequence.

For example:
- `[]` specifies no modifications at all
- `[with_changed_bound("CO2t", lb=0, ub=10)]` limits the CO2 transport
- `[with_changed_bound("CO2t", lb=0, ub=2), with_changed_bound("O2t", lb=0, ub=100)]`
  severely limits the CO2 transport _and_ slightly restricts the transport of
  O2

!!! note "Variants are single-parameter model-transforming functions"
	Because the variants are just generators of single parameter functions
	that take the model and return its modified version, you can also use
	`identity` to specify a variant that does nothing -- `[identity]` is
	perfectly same as `[]`

The shape of the variants array is important too, because it is precisely
retained in the result (just as with `pmap`). If you pass in a matrix of
variants, you will receive a matrix of analysis results of the same size. That
can be exploited for easily exploring many combinations of possible model
properties. Let's try exploring a "cube" of possible restricted reactions:

```julia
using IterTools # for cartesian products

res = screen_variants(m,
    [
        # for each variant we restricts 2 reactions
        [with_changed_bound(r1, lb=-3, ub=3), with_changed_bound(r2, lb=-1, ub=1)]

        # the reaction pair will be chosen from a cartesian product
        for (r1,r2) in product(
            ["H2Ot", "CO2t", "O2t", "NH4t"], # of this set of transport reactions
            ["EX_h2o_e", "EX_co2_e", "EX_o2_e", "EX_nh4_e"], # and this set of exchanges
        )
    ],
    m -> flux_balance_analysis_dict(m, Tulip.Optimizer)["BIOMASS_Ecoli_core_w_GAM"],
)
```

As a result, we will receive a full matrix of the biomass productions:
```
4×4 Matrix{Float64}:
 0.407666  0.454097  0.240106  0.183392
 0.407666  0.485204  0.24766   0.183392
 0.314923  0.319654  0.24766   0.183392
 0.407666  0.485204  0.24766   0.183392
```
Notably, this shows that O2 transport and NH4 exchange may be serious
bottlenecks for biomass production.

For clarity, you may always annotate the result by zipping it with the
specification structure you have used and collecting the data:
``` julia
collect(zip(
    product(
        ["H2Ot", "CO2t", "O2t", "NH4t"],
        ["EX_h2o_e", "EX_co2_e", "EX_o2_e", "EX_nh4_e"],
    ),
    res,
))
```
...which gives the following annotated result:
```
4×4 Matrix{Tuple{Tuple{String, String}, Float64}}:
 (("H2Ot", "EX_h2o_e"), 0.407666)  (("H2Ot", "EX_co2_e"), 0.454097)  (("H2Ot", "EX_o2_e"), 0.240106)  (("H2Ot", "EX_nh4_e"), 0.183392)
 (("CO2t", "EX_h2o_e"), 0.407666)  (("CO2t", "EX_co2_e"), 0.485204)  (("CO2t", "EX_o2_e"), 0.24766)   (("CO2t", "EX_nh4_e"), 0.183392)
 (("O2t", "EX_h2o_e"), 0.314923)   (("O2t", "EX_co2_e"), 0.319654)   (("O2t", "EX_o2_e"), 0.24766)    (("O2t", "EX_nh4_e"), 0.183392)
 (("NH4t", "EX_h2o_e"), 0.407666)  (("NH4t", "EX_co2_e"), 0.485204)  (("NH4t", "EX_o2_e"), 0.24766)   (("NH4t", "EX_nh4_e"), 0.183392)
```

This may be easily used for e.g. scrutinizing all possible reaction pairs, to
find the ones that are redundant and not.

There are many other variant "specifications" to choose from. You may use
[`with_added_reactions`](@ref), [`with_removed_reactions`](@ref),
[`with_removed_metabolites`](@ref), and others. Function reference contains a
complete list; as a convention, names of the specifications all start with
`with_`.

## Writing custom variant functions

It is actually very easy to create custom specifications that do any
modification that you can implement, to be later used with
[`screen_variants`](@ref) and [`screen`](@ref).

Generally, the "specifications" are supposed to return a _function_ that
creates a modified copy of the model. The copy of the model may be shallow, but
the functions should always prevent modifying the original model structure --
`screen` is keeping a single copy of the original model at each worker to
prevent unnecessary bulk data transport, and if that is changed in-place, all
following analyses of the model will work on inconsistent data, usually
returning wrong results (even randomly changing ones, because of the
asynchronous nature of [`screen`](@ref) execution).

Despite of that, writing a modification is easy. The simplest modification that
"does nothing" (isomorphic to standard `identity`) can be formatted as follows:

```julia
with_no_change = model -> model
```

The modifications may change the model, provided it is copied properly. The
following modification will remove a reaction called "O2t", effectively
removing the possibility to transport oxygen. We require a specific type of
model where this change is easy to perform (generally, not all variants may be
feasible on all model types).

```julia
with_disabled_oxygen_transport = (model::ObjectModel) -> begin

    # make "as shallow as possible" copy of the `model`.
    # Utilizing `deepcopy` is also possible, but inefficient.
    new_model = copy(model)
    new_model.reactions = copy(model.reactions)

    # remove the O2 transport from the model copy
    delete!(new_model.reactions, "O2t")

    return new_model #return the newly created variant
end
```

Finally, the whole definition may be parameterized as a normal function. The
following variant removes any user-selected reaction:

```julia
with_disabled_reaction(reaction_id) = (model::ObjectModel) -> begin
    new_model = copy(model)
    new_model.reactions = copy(model.reactions)
    delete!(new_model.reactions, reaction_id) # use the parameter from the specification
    return new_model
end
```

In turn, these variants can be used in [`screen_variants`](@ref) just as we
used [`with_changed_bound`](@ref) above:

```julia
screen_variants(
    m,    # the model for screening
    [
        [with_no_change],
        [with_disabled_oxygen_transport],
        [with_disabled_reaction("NH4t")],
    ],
    m -> flux_balance_analysis_dict(m, Tulip.Optimizer)["BIOMASS_Ecoli_core_w_GAM"],
)
```

That should get you the results for all new variants of the model:
```
3-element Vector{Float64}:
 0.8739215022674809
 0.21166294865468896
 1.2907224478973395e-15
```

!!! warning "Custom variants with distributed processing"
    If using distributed evaluation, remember the variant-generating functions
    need to be defined on all used workers (generating the variants in parallel
    on the workers allows COBREXA to run the screening process very
    efficiently, without unnecessary sending of bulk model data). Prefixing the
    definition with `@everywhere` is usually sufficient for that purpose.

## Passing extra arguments to the analysis function

Some analysis functions may take additional arguments, which you might want to
vary for the analysis. `modifications` argument of
[`flux_balance_analysis_dict`](@ref) is one example of such argument, allowing
you to specify details of the optimization procedure.

[`screen`](@ref) function allows you to do precisely that -- apart from
`variants`, you may also specify an array of `args` of the same shape as
`variants`, the entries of which will get passed together with the generated
model variants to your specified analysis function. If either of the arguments
is missing (or set to `nothing`), it is defaulted to "no modifications" or "no
arguments".

The arguments _must_ be tuples; you may need to make 1-tuples from your data
(e.g. using `(value,)`) if you want to pass just a single argument.

Let's try to use that functionality for trying to find a sufficient amount of
iterations needed for Tulip solver to find a feasible solution:

```julia
screen(m,
    args = [(i,) for i in 5:15],  # the iteration counts, packed in 1-tuples
    analysis = (m,a) -> # `args` elements get passed as the extra parameter here
        flux_balance_analysis_vec(m,
            Tulip.Optimizer;
            modifications=[change_optimizer_attribute("IPM_IterationsLimit", a)],
        ),
)
```

From the result, we can see that Tulip would need at least 14 iterations to
find a feasible region:

```
11-element Vector{Union{Nothing, Vector{Float64}}}:
 nothing
 nothing
 nothing
 nothing
 nothing
 nothing
 nothing
 nothing
 nothing
 [7.47738193404817, 1.8840414375838503e-8, 4.860861010127701, -16.023526104614593, … ]
 [7.47738193404817, 1.8840414375838503e-8, 4.860861010127701, -16.023526104614593, … ]
```
