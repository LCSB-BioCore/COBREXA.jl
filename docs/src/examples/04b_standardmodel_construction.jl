# # Model construction and modification

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats; and convert
# these into `ObjectModel`s. However, it is also possible to construct models
# from scratch, and modify existing models. This will be demonstrated
# here.

using COBREXA

# In `COBREXA`, model construction is primarily supported through `ObjectModel`s.
# To begin, create an empty `ObjectModel`.

model = ObjectModel()

# Next, genes, metabolites and reactions need to be added to the model.

# ### Add genes to the model
gene_list = [Gene(string("g", num)) for num = 1:8]

#md # !!! warning "Warning: Don't accidentally overwrite the generic accessors"
#md #       It may be tempting to call a variable `genes`, `metabolites`, or
#md #       `variables`. However, these names conflict with generic accessors
#md #       functions and will create problems downstream.

add_genes!(model, gene_list)

# ### Add metabolites to the model

metabolite_list = [Metabolite(string("m", num)) for num = 1:4]

metabolite_list[1].formula = "C6H12O6" # can edit metabolites, etc. directly

add_metabolites!(model, metabolite_list)

# ### Add reactions to the model

r_m1 = ReactionBidirectional("EX_m1", Dict("m1" => -1.0)) # exchange reaction: m1 <-> (is the same as m1 ↔ nothing)
r1 = ReactionForward("r1", Dict("m1" => -1.0, "m2" => 1.0))
r1.gene_associations = [Isozyme(["g1", "g2"]), Isozyme(["g3"])] # add some gene reaction rules
r2 = ReactionBackward("r2", Dict("m2" => -1.0, "m1" => 1.0))
r3 = ReactionBidirectional("r3", Dict("m2" => -1.0, "m3" => 1.0))
r4 = ReactionForward("r3", Dict("m2" => -1.0, "m4" => 1.0))
r_m3 = ReactionBidirectional("r3", Dict("m3" => -1.0))
r_m4 = ReactionForward("r3", Dict("m4" => -1.0))
r5 = ReactionForward("r5", Dict("m4" => -1.0, "m2" => 1.0))

add_reactions!(model, [r1, r2, r3, r_m1, r4, r_m3, r_m4, r5])

m1 = metabolite_list[1]
m2 = metabolite_list[2]
m3 = metabolite_list[3]
m4 = metabolite_list[4]

model.reactions["r4"].gene_associations =
    [Isozyme(x) for x in [["g5"], ["g6", "g7"], ["g8"]]]

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

# Likewise, reactions and metabolites can also be deleted.

remove_metabolite!(model, "m1")
model
