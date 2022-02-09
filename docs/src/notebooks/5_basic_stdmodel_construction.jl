# # Model construction and modification

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats; and convert
# these into `StandardModel`s. However, it is also possible to construct models
# from scratch, and modify existing models. This will be demonstrated
# here.

# ## Model construction

# In `COBREXA`, model construction is primarily supported through `StandardModel`s.
# To begin, create an empty `StandardModel`.

using COBREXA

model = StandardModel("FirstModel") # assign model id = "FirstModel"

# Next, genes, metabolites and reactions need to be added to the model.

# ### Add genes to the model
gene_list = [Gene(string("g", num)) for num = 1:8]

#md # !!! warning "Warning: Don't accidentally overwrite the generic accessors"
#md #       It may be tempting to call a variable `genes`, `metabolites`, or
#md #       `reactions`. However, these names conflict with generic accessors
#md #       functions and will create problems downstream.

add_genes!(model, gene_list)

# ### Add metabolites to the model
metabolite_list = [Metabolite(string("m", num)) for num = 1:4]

metabolite_list[1].formula = "C6H12O6" # can edit metabolites, etc. directly

add_metabolites!(model, metabolite_list)

# ### Add reactions to the model

# There are two ways to create and add reactions to a model.
# These are using functions, or macros.

r_m1 = Reaction("EX_m1", Dict("m1" => -1.0), :bidirectional) # exchange reaction: m1 <-> (is the same as m1 ↔ nothing)
r1 = Reaction("r1", Dict("m1" => -1.0, "m2" => 1.0), :forward)
r1.grr = [["g1", "g2"], ["g3"]] # add some gene reaction rules
r2 = Reaction("r2", Dict("m2" => -1.0, "m1" => 1.0), :reverse)
r3 = Reaction("r3", Dict("m2" => -1.0, "m3" => 1.0), :bidirectional)

add_reactions!(model, [r1, r2, r3, r_m1]) # function approach

m1 = metabolite_list[1]
m2 = metabolite_list[2]
m3 = metabolite_list[3]
m4 = metabolite_list[4]

@add_reactions! model begin # macro approach
    "r4", m2 → m4, 0, 1000
    "r_m3", m3 ↔ nothing, -1000, 1000
    "r_m4", m4 → nothing
    "r5", m4 → m2
end

model.reactions["r4"].grr = [["g5"], ["g6", "g7"], ["g8"]]

#md # !!! note "Note: Writing unicode arrows"
#md #     The reaction arrows can be easily written by using the `LaTeX`
#md #     completions built into Julia shell (and many Julia-compatible
#md #     editors). You can type:
#md #
#md #     - `→` as `\rightarrow` (press `Tab` to complete)
#md #     - `←` as `\leftarrow`
#md #     - `↔` as `\leftrightarrow`

# The constructed model can now be inspected.
model

# ## Modifying existing models

# It is also possible to modify a model by deleting certain genes.
# This is simply achieved by calling `remove_genes!`.

remove_genes!(model, ["g1", "g2"]; knockout_reactions = false)
model
