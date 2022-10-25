
# Extending the models

To simplify doing (and undoing) simple modifications to the existing model
structure, COBREXA.jl supports a class of model _wrappers_, which are basically
small layers that add or change the functionality of a given base models.

Types [`Serialized`](@ref), [`CoreCoupling`](@ref), [`SMomentModel`](@ref), and
[`GeckoModel`](@ref) all work in this manner -- add some extra functionality to
the "base". Technically, they are all subtypes of the abstract type
[`ModelWrapper`](@ref), which itself is a subtype of [`AbstractMetabolicModel`](@ref)
and can thus be used in all standard analysis functions.  Similarly, the model
wraps can be stacked -- it is easy to e.g. serialize a [`GeckoModel`](@ref), or
to add coupling to an existing [`SMomentModel`](@ref).

As the main benefit of the approach, creating model variants using the wrapper
approach is usually more efficient than recomputing the models in place. The
wrappers are thin, and if all values can get computed and materialized only once
the model data is actually needed, we may save a great amount of computing
power.

At the same time, since the original model stays unchanged (and may even be
immutable), undoing the modifications caused by the wrapper is extremely easy
and fast -- we just discard the wrapper.

## Writing a model wrapper

Creating a model wrapper structure is simple -- by declaring it a subtype of
[`ModelWrapper`](@ref) and implementing a single function
[`unwrap_model`](@ref), we get default implementations of all accessors that
should work for any [`AbstractMetabolicModel`](@ref).

As a technical example, we may make a minimal model wrapper that does not do
anything:

```julia
struct IdentityWrap <: ModelWrapper
    mdl::AbstractMetabolicModel
end

COBREXA.unwrap_model(x::IdentityWrap) = x.mdl
```

This is instantly usable in all analysis functions, although there is no
actual "new" functionality:

```julia
m = IdentityWrap(load_model("e_coli_core.xml"))
flux_balance_analysis_vec(m, GLPK.Optimizer)
```

To modify the functionality, we simply add specific methods for accessors that
we want modified, such as [`bounds`](@ref), [`stoichiometry`](@ref) and
[`objective`](@ref). We demonstrate that on several examples below.

## Example 1: Slower model

Here, we construct a type `RateChangedModel` that has all bounds multiplied by
a constant factor. This can be used to e.g. simulate higher or lower abundance
of certain organism in a model.

```julia
struct RateChangedModel <: ModelWrapper
    factor::Float64
    mdl::AbstractMetabolicModel
end
```

The overloaded accessors typically reach for basic information into the "inner"
wrapped model, and modify them in a certain way.

```julia
COBREXA.unwrap_model(x::RateChangedModel) = x.mdl
function COBREXA.bounds(x::RateChangedModel)
    (l, u) = bounds(x.mdl) # extract the original bounds
    return (l .* x.factor, u .* x.factor) # return customized bounds
end
```

To make a 2 times faster or slower model from a base model, we can run:
```julia
faster_e_coli = RateChangedModel(2.0, load_model("e_coli_core.xml"))
slower_e_coli = RateChangedModel(1/2, load_model("e_coli_core.xml"))
```

## Example 2: Leaky model

As the second example, we construct a hypothetical model that is "leaking" all
metabolites at once at a constant fixed small rate. Again, the modification is
not quite realistic, but may be useful to validate the mathematical robustness
of the models.

```julia
struct LeakyModel <: ModelWrapper
    leaking_metabolites::Vector{String}
    leak_rate::Float64
    mdl::AbstractMetabolicModel
end
```

Technically, we implement the leaks by adding an extra reaction bounded to the
precise `leak_rate`, which permanently removes all metabolites. That is done by
modifying the reaction list, stoichiometry, and bounds:

```julia
COBREXA.unwrap_model(x::LeakyModel) = x.mdl
COBREXA.n_reactions(x::LeakyModel) = n_reactions(x.mdl) + 1
COBREXA.reactions(x::LeakyModel) = [reactions(x.mdl); "The Leak"]
COBREXA.stoichiometry(x::LeakyModel) = [stoichiometry(x.mdl) [m in x.leaking_metabolites ? -1.0 : 0.0 for m = metabolites(x.mdl)]]
function COBREXA.bounds(x::LeakyModel)
    (l, u) = bounds(x.mdl)
    return ([l; x.leak_rate], [u; x.leak_rate])
end
```

To make the wrapper complete and consistent, we also have to modify the
accessors that depend on correct sizes of the model items.

```julia
COBREXA.objective(x::LeakyModel) = [objective(x.mdl); 0]
COBREXA.reaction_flux(x::LeakyModel) = [reaction_flux(x.mdl); zeros(1, n_reactions(x.mdl))]
COBREXA.coupling(x::LeakyModel) = [coupling(x.mdl) zeros(n_coupling_constraints(x.mdl))]
```
(Among other, we modified the [`reaction_flux`](@ref) so that all analysis
methods ignore the leak reaction.)

Now, any model can be made to lose some chosen metabolites as follows:
```julia
leaks = ["M_o2_c", "M_pi_c", "M_glx_c"]
leaky_e_coli = LeakyModel(leaks, 5, load_model("e_coli_core.xml"))
```

## Example 3: Combining the wrappers

With both wrappers implemented individually, it is easy to combine them by
re-wrapping. We can easily create a model that is slowed down and moreover
leaks the metabolites as follows:
```julia
leaky_slow_e_coli = LeakyModel(leaks, 5, RateChangedModel(1/2, load_model("e_coli_core.xml")))
```

As with all wrapping operations, take care about the exact order of applying
the wraps. The other combination of the model wraps differs by also changing
the rate of the metabolite leaks, which did not happen with the
`leaky_slow_e_coli` above:
```julia
slowly_leaking_slow_e_coli = RateChangedModel(1/2, LeakyModel(leaks, 5, load_model("e_coli_core.xml")))
```

Expectably, the model can be solved with standard functions:
```julia
v = flux_balance_analysis_dict(slowly_leaking_slow_e_coli, GLPK.Optimizer)
v["R_BIOMASS_Ecoli_core_w_GAM"]  # prints out ~0.38
```
