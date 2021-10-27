# # Basic usage of `StandardModel`

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this tutorial we will use `COBREXA`'s `StandardModel` and functions that
# specifically operate on it. As usual we will use the toy model of *E. coli*
# for demonstration.

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

using COBREXA

# ## Loading a model

model = load_model(StandardModel, "e_coli_core.json") # we specifically want to load a StandardModel from the model file

#md # !!! note "Note: Loading `StandardModel`s implicitly uses `convert`"
#md #       When using `load_model(StandardModel, file_location)` the model at
#md #       `file_location` is first loaded into its inferred format and is then
#md #       converted to a `StandardModel` using the generic accessor interface.
#md #       Thus, data loss may occur. Always check your model to ensure that
#md #       nothing important has been lost.

#nb # When using `load_model(StandardModel, file_location)` the model at
#nb # `file_location` is first loaded into its inferred format and is then
#nb # converted to a `StandardModel` using the generic accessor interface.
#nb # Thus, data loss may occur. Always check your model to ensure that
#nb # nothing important has been lost.

# ## Basic analysis

# As before, for optimization based analysis we need to load an optimizer. Here we
# will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) to solve the linear
# programs of this tutorial. Refer to the basic constraint-based analysis
# tutorial for more informaiton.

# All the normal analysis functions work on `StandardModel`, due to it also
# having the same generic accessor interface as all the other model types.

using Tulip

fluxes = flux_balance_analysis_dict(
    model,
    Tulip.Optimizer;
    modifications = [
        change_objective("BIOMASS_Ecoli_core_w_GAM"),
        change_constraint("EX_glc__D_e"; lb = -12, ub = -12),
        change_constraint("EX_o2_e"; lb = 0, ub = 0),
    ],
)

# This is not very exciting yet, since every other model type can also do this.
# However, deeper inspection of flux results is possible when using
# `StandardModel`.

# ## Inspecting the flux solution: `atom_fluxes`

# It is sometimes interesting to keep track of the atoms entering and leaving
# the system through boundary reactions. This can be inspected by calling
# [`atom_fluxes`](@ref). That gives you the flux of individual atoms getting
# consumed and produced by all reactions, based on `fluxes`. We erase the
# reaction that consumes the atoms for creating biomass, to see how much mass
# the "rest" of the reaction produces for it:

fluxes_without_biomass = copy(fluxes);
delete!(fluxes_without_biomass, "BIOMASS_Ecoli_core_w_GAM");
atom_fluxes(model, fluxes_without_biomass)

# ## Inspecting the flux solution: `metabolite_fluxes`

# Another useful flux result analysis function is [`metabolite_fluxes`](@ref).
# This function gives an overview of reactions consuming and producing each
# metabolite.

consuming, producing = metabolite_fluxes(model, fluxes)

consuming["atp_c"] # reactions consuming `atp_c`

# ## Internals of `StandardModel`

# Another benefit of `StandardModel` is that it supports a richer internal
# infrastructure that can be used to manipulate internal model attributes in a
# systematic way. Specifically, the genes, reactions, and metabolites with of a
# model each have a type. This is particularly useful when modifying or even
# constructing a model from scratch.

# ## `Gene`s, `Reaction`s, and `Metabolite`s

# `StandardModel` is composed of ordered dictionaries of `Gene`s, `Metabolite`s
# and `Reaction`s. Ordered dictionaries are used because the order of the
# reactions and metabolites are important for constructing a stoichiometric
# matrix since the rows and columns should correspond to the order of the metabolites
# and reactions returned by calling the accessors `metabolites` and `reactions`.

# Each `StandardModel` is composed of the following fields:

fieldnames(StandardModel) # fields of a StandardModel

# The `:genes` field of a `StandardModel` contains an ordered dictionary of gene ids mapped to `Gene`s.

model.genes # the keys of this dictionary are the same as genes(model)

# The `Gene` type is a struct that can be used to store information about genes
# in a `StandardModel`. Each `Gene` is composed of the following fields:

fieldnames(Gene)

#md # !!! tip "Tip: Use <tab> complete to explore the structure of types"
#md #       Use <tab> to quickly explore the fields of a struct. For example,
#md #       Gene.<tab> will list all the fields shown above.

#nb # Use <tab> to quickly explore the fields of a struct. For example,
#nb # Gene.<tab> will list all the fields shown above.

# The keys used in the ordered dictionaries in
# `model.genes` are the ids returned using the generic accessor `genes`. `Gene`s
# have pretty printing, as demonstrated below for a random gene drawn from the
# model:

random_gene_id = genes(model)[rand(1:n_genes(model))]
model.genes[random_gene_id]

# The same idea holds for both metabolites (stored as `Metabolite`s) and
# reactions (stored as `Reaction`s). This is demonstrated below.

random_metabolite_id = metabolites(model)[rand(1:n_metabolites(model))]
model.metabolites[random_metabolite_id]
#
random_reaction_id = reactions(model)[rand(1:n_reactions(model))]
model.reactions[random_reaction_id]

# `StandardModel` can be used to build your own metabolic model or modify an
# existing one. One of the main use cases for `StandardModel` is that it can be
# used to merge multiple models or parts of multiple models together. Since the
# internals are uniform inside each `StandardModel`, attributes of other model
# types are squashed into the required format (using the generic accessors).
# This ensures that the internals of all `StandardModel`s are the same -
# allowing easy systematic evaluation.

#md # !!! warning "Warning: Combining models with different namespaces is tricky"
#md #       Combining models that use different namespaces requires care.
#md #       For example, in some models the water exchange reaction is called
#md #       `EX_h2o_e`, while in others it is called `R_EX_h2o_s`. This needs to
#md #       manually addressed (for now) to prevent duplicate, e.g. reactions,
#md #       from being added.

# ## Checking the internals of `StandardModel`s: `annotation_index`

# Often when models are automatically reconstructed duplicate genes, reactions
# or metabolites end up in a model. `COBREXA` exports `annotation_index` to
# check for cases where the id of a struct may be different, but the annotations
# the same (possibly suggesting a duplication). `annotation_index` builds a
# dictionary mapping annotation features to the ids of whatever struct you are
# inspecting. This makes it easy to find structs that share certain annotation features.

rxn_annotations = annotation_index(model.reactions)
#
rxn_annotations["ec-code"]

# The `annotation_index` function can also be used on `Reaction`s and
# `Gene`s in the same way.

# ## Checking the internals of `StandardModel`s: `check_duplicate_reaction`

# Another useful function is `check_duplicate_reaction`, which checks for
# reactions that have duplicate (or similar) reaction equations.

pgm_duplicate = Reaction()
pgm_duplicate.id = "pgm2" # Phosphoglycerate mutase
pgm_duplicate.metabolites = Dict{String,Float64}("3pg_c" => 1, "2pg_c" => -1)
pgm_duplicate
#
check_duplicate_reaction(pgm_duplicate, model.reactions; only_metabolites = false) # can also just check if only the metabolites are the same but different stoichiometry is used

# ## Checking the internals of `StandardModel`s: `reaction_mass_balanced`

# Finally, [`reaction_mass_balanced`](@ref) can be used to check if a reaction is mass
# balanced based on the formulas of the reaction equation.

rxn_dict = Dict{String,Float64}("3pg_c" => 1, "2pg_c" => -1, "h2o_c" => 1)
reaction_mass_balanced(model, rxn_dict)

# Now to determine which atoms are unbalanced, you can use `reaction_atom_balance`
reaction_atom_balance(model, rxn_dict)

# Note, since `pgm_duplicate` is not in the model, we cannot use the other variants of this
# function because they find the reaction equation stored inside the `model`.
