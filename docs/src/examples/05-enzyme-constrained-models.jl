# # Enzyme constrained models 

using COBREXA

# Here we will construct an enzyme constrained variant of the *E. coli* "core"
# model. We will need the model, which we can download if it is not already present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# Enzyme constrained models require parameters that are usually not used by
# conventional constraint based models. These include reaction specific turnover
# numbers, molar masses of enzymes, and capacity bounds. 

import AbstractFBCModels as A

### Reaction turnover numbers

# Here we will use randomly generated turnover numbers for simplicity. Each
# reaction in a constraint-based model usually has gene reaction rules
# associated with it. These typically take the form of, possibly multiple,
# isozymes that can catalyze a reaction. A turnover number needs to be assigned
# to each isozyme, as shown below.

using Random

reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue #src
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,X.Isozyme}())
        # each isozyme gets a unique name
        d["isozyme_"*string(i)] = X.Isozyme( # Isozyme struct is defined by COBREXA 
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = 100 + 50 * rand(), # forward reaction turnover number
            kcat_backward = 100 + 50 * rand(), # reverse reaction turnover number
        )
        d["isozyme_"*string(i)] = X.Isozyme( #src
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), #src
            kcat_forward = ecoli_core_reaction_kcats[rid][i][1], #src
            kcat_backward = ecoli_core_reaction_kcats[rid][i][2], #src
        )
    end
end

#!!! warning "Turnover number units"
#    Take care with the units of the turnover numbers. In literature they are 
#    usually reported in 1/s. However, flux units are typically mmol/gDW/h, 
#    suggesting that you should rescale the turnover numbers to 1/h if you 
#    want to use the traditional flux units. 

### Enzyme molar masses

# Similarly, we will randomly generate enzyme molar masses for use in the enzyme
# constrained model.

gene_molar_masses = Dict(gid => 100 + 40 * rand() for gid in A.genes(model))
gene_molar_masses = ecoli_core_gene_product_masses #src

#!!! warning "Molar mass units"
#    Take care with the units of the molar masses. In literature they are 
#    usually reported in Da or kDa (g/mol). However, as noted above, flux 
#    units are typically mmol/gDW/h. Since the enzyme kinetic equation is 
#    `v = k * e`, where `k` is the turnover number, it suggests that the 
#    enzyme variable will have units of mmol/gDW. The molar masses come 
#    into play when setting the capacity limitations, e.g. usually a sum 
#    over all enzymes weighted by their molar masses: `e * mm`. Thus, if 
#    your capacity limitation has units of g/gDW, then the molar masses 
#    must have units of g/mmol.  

### Capacity limitation

# The capacity limitation usually denotes an upper bound of protein available to
# the cell. 

total_enzyme_capacity = 100.0 # g enzyme/gDW

### Running a basic enzyme constrained model

m = build_enzyme_constrained_model(
    model,
    reaction_isozymes,
    gene_molar_masses,
    [("total_proteome_bound", A.genes(model), total_enzyme_capacity)],
)

m.fluxes.EX_glc__D_e.bound = (-1000.0, 1000.0) # undo glucose important bound
m.fluxes.GLCpts.bound = (-1.0, 12.0) # undo glucose important bound
for (k, v) in m.enzymes
    if k == :b2779
        v.bound = (0.01, 0.06)
    else
        v.bound = (0.0, 1.0)
    end
end

sol = optimized_constraints(
    m;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications =[set_optimizer_attribute("IPM_IterationsLimit", 1000)]
)


### Building the model incrementally
