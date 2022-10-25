
# Working with custom models

It may happen that the intuitive representation of your data does not really
match what is supported by a given COBRA package. COBREXA.jl attempts to avoid
this problem by providing a flexible framework for containing any data
structure that can, somehow, represent the constraint-based model.

The task of having such a polymorphic model definition can be split into 2
separate concerns:

- How to allow the analysis functions to gather the required information from
  any user-specified model data structure?
- How to make the reconstruction functions (i.e., reaction or gene deletions)
  work properly on any data structure?

To solve the first concern, COBREXA.jl specifies a set of generic accessors
that work over the abstract type [`AbstractMetabolicModel`](@ref). To use your data
structure in a model, you just make it a subtype of [`AbstractMetabolicModel`](@ref)
and overload the required accessors. The accessors are functions that extract
some relevant information, such as [`stoichiometry`](@ref) and
[`bounds`](@ref), returning a fixed simple data type that can be further used
by COBREXA.  You may see a complete list of accessors
[here](../functions.md#Base-Types).

A good solution to the second concern is a slightly more involved, as writing
generic data modifiers is notoriously hard. Still, there is support for easily
making small changes to the model using the modifications system, with
functions such as [`with_added_reactions`](@ref) and
[`with_changed_bound`](@ref).

## Writing the generic accessors

Let's write a data structure that represents a very small model that contains N
metabolites that are converted in a circle through N linear, coupled reactions.
(E.g., for N=3, we would have a conversion of metabolites A, B and C ordered as
A → B → C → A.) This may be useful for testing purposes; we will use it for a
simple demonstration.

The whole model can thus be specified with a single integer N that represents
the length of the reaction cycle:

```julia
struct CircularModel <: AbstractMetabolicModel
    size::Int
end
```

First, define the reactions and metabolites:

```julia
COBREXA.n_reactions(m::CircularModel) = m.size
COBREXA.n_metabolites(m::CircularModel) = m.size

COBREXA.reactions(m::CircularModel) = ["rxn$i" for i in 1:n_reactions(m)]
COBREXA.metabolites(m::CircularModel) = ["met$i" for i in 1:n_metabolites(m)]
```

It is useful to re-use the already defined functions, as that improves the code
maintainability.

We can continue with the actual linear model properties:

```julia
function COBREXA.objective(m::CircularModel)
    c = spzeros(n_reactions(m))
    c[1] = 1 #optimize the first reaction
    return c
end

COBREXA.bounds(m::CircularModel) = (
    zeros(n_reactions(m)), # lower bounds
    ones(n_reactions(m)), # upper bounds
)

function COBREXA.stoichiometry(m::CircularModel)
    nr = n_reactions(m)
    stoi(i,j) =
        i == j ? 1.0 :
        (i % nr + 1) == j  ? -1.0 :
        0.0

    sparse([stoi(i,j) for i in 1:nr, j in 1:nr])
end
```

You may check that the result now works just as with [`CoreModel`](@ref) and
[`StandardModel`](@ref):

```julia
julia> m = CircularModel(5)
Metabolic model of type CircularModel

  1.0  -1.0    ⋅     ⋅     ⋅
   ⋅    1.0  -1.0    ⋅     ⋅
   ⋅     ⋅    1.0  -1.0    ⋅
   ⋅     ⋅     ⋅    1.0  -1.0
 -1.0    ⋅     ⋅     ⋅    1.0
Number of reactions: 5
Number of metabolites: 5

```

This interface is sufficient to run most of the basic analyses, especially the flux balance finding ones:

```julia
julia> flux_balance_analysis_dict(m, Tulip.Optimizer)
Dict{String, Float64} with 5 entries:
  "rxn5" => 1.0
  "rxn2" => 1.0
  "rxn1" => 1.0
  "rxn3" => 1.0
  "rxn4" => 1.0

```

## Writing generic model modifications

The custom model structure can also be made compatible with many of the
existing variant-generating functions and analysis modifiers.

The functions prepared for use as "variants" in [`screen`](@ref), usually
prefixed by `with_`, have their generic variants that only call simpler,
overloadable functions for each specific model. This choice is based on the
overloading dispatch system of Julia. For
example,[`with_removed_metabolites`](@ref) is implemented very generically by
reducing the problem to some specific [`remove_metabolites`](@ref) functions
selected by the dispatch, as follows:

```julia
with_removed_metabolites(args...; kwargs...) =
    m -> remove_metabolites(m, args...; kwargs...)
```

To be able to use [`with_removed_metabolites`](@ref) in your model, we can just
overload the actual inner function. For the simple circular model, the
modification might as well look like this:

```julia
COBREXA.remove_metabolites(m::CircularModel, n::Int) =
    return CircularModel(m.size - n)
```

!!! danger "Functions that generate model variants must be pure"
    Notice that the function is "pure", i.e., does not make any in-place
    modifications to the original model structure. That property is required
    for [`screen`](@ref) and other functions to properly and predictably apply
    the modifications to the model. To expose potential in-place modifications
    to your codebase, you should instead overload the "bang" counterpart of
    remove metabolites, called [`remove_metabolites!`](@ref).
