# # Model construction and modification

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats; and convert
# these into `StandardModel`s. However, it is also possible to construct models
# from scratch, join models together, and modify exising models. This will be demonstrated
# here.

# ## Model construction 

# In `COBREXA`, model construction is primarily supported through `StandardModel`s.  
# To begin, create an empty `StandardModel`.

using COBREXA

model = StandardModel("FirstModel") # assign model id = "FirstModel"

# Next, genes, metabolites and reactions need to be added to the model.

# ### Add genes to the model
gene_list = [Gene(string("g",num)) for num = 1:8]

#md # !!! warning "Warning: Don't accidentally overwrite the generic accessors"
#md #       It may be tempting to call a variable `genes`, `metabolites`, or 
#md #       `reactions`. However, these names conflict with generic accessors
#md #       functions and will create problems downstream.

add!(model, gene_list)

# ### Add metabolites to the model
metabolite_list = [Metabolite(string("m", num)) for num = 1:4]

metabolite_list[1].formula = "C6H12O6" # can edit metabolites, etc. directly

add!(model, metabolite_list)

# ### Add reactions to the model

# There are two ways to create reactions. Using a function and using a macro.

r_m1 = Reaction("EX_m1", Dict("m1" => -1.0), :bidirectional) # exchange reaction: m1 <-> ∅ (nothing)
r1 = Reaction("r1", Dict("m1" => -1.0, "m2" => 1.0), :forward)
r2 = Reaction("r2", Dict("m2" => -1.0, "m1" => 1.0), :backward)
r3 = Reaction("r3", Dict("m2" => -1.0, "m3" => 1.0), :bidirectional)

add!(model, [r1, r2, r3, r_m1]) 

m1 = metabolite_list[1]
m2 = metabolite_list[2]
m3 = metabolite_list[3]
m4 = metabolite_list[4]

@add_reactions! model begin
    r4, m2 ⟶ m4, 0, 1000
    r_m3, m3 ⟷ ∅, -1000, 1000
    r_m4, m4 ⟶ ∅
end
