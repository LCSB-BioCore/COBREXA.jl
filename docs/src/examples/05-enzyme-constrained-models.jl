
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Enzyme constrained models

using COBREXA

# Here we will construct an enzyme constrained variant of the *E. coli* "core"
# model. We will need the model, which we can download if it is not already present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use Tulip here:

import AbstractFBCModels as A
import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# Enzyme constrained models require parameters that are usually not used by
# conventional constraint based models. These include reaction specific turnover
# numbers, molar masses of enzymes, and capacity bounds.

# ## Reaction turnover numbers

# Enzyme constrained models require reaction turnover numbers, which are often
# isozyme specfic. Many machine learning tools, or experimental data sets, can
# be used to estimate these parameters.

# ```@raw html
# <details><summary>Reaction turnover numbers</summary>
# ??? how to load this?
# </details>
# ```

# Each reaction in a constraint-based model usually has gene reaction rules
# associated with it. These typically take the form of, possibly multiple,
# isozymes that can catalyze a reaction. A turnover number needs to be assigned
# to each isozyme, as shown below.

reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
        # each isozyme gets a unique name
        d["isozyme_"*string(i)] = Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = ecoli_core_reaction_kcats[rid] * 3.6, # forward reaction turnover number units = 1/h
            kcat_reverse = ecoli_core_reaction_kcats[rid] * 3.6, # reverse reaction turnover number units = 1/h
        )
    end
end

#!!! warning "Turnover number units"
#    Take care with the units of the turnover numbers. In literature they are
#    usually reported in 1/s. However, flux units are typically mmol/gDW/h,
#    suggesting that you should rescale the turnover numbers to 1/h if you
#    want to use the conventional flux units.

### Enzyme molar masses

# We also require the mass of each enzyme, to properly weight the contribution
# of each flux/isozyme in the capacity bound(s). These data can typically be
# found in uniprot.

# ```@raw html
# <details><summary>Reaction turnover numbers</summary>
# ??? how to load this?
# </details>
# ```

gene_product_molar_masses = ecoli_core_gene_product_masses

#!!! warning "Molar mass units"
#    Take care with the units of the molar masses. In literature they are
#    usually reported in Da or kDa (g/mol). However, as noted above, flux
#    units are typically mmol/gDW/h. Since the enzyme kinetic equation is
#    `v = k * e`, where `k` is the turnover number, it suggests that the
#    enzyme variable will have units of mmol/gDW. The molar masses come
#    into play when setting the capacity limitations, e.g. usually a sum
#    over all enzymes weighted by their molar masses: `e * mm`. Thus, if
#    your capacity limitation has units of g/gDW, then the molar masses
#    must have units of g/mmol (= kDa).

### Capacity limitation

# The capacity limitation usually denotes an upper bound of protein available to
# the cell.

total_enzyme_capacity = 50.0 # mg enzyme/gDW

### Running a basic enzyme constrained model

# With all the parameters specified, we can directly use the enzyme constrained
# convenience function to run enzyme constrained FBA in one shot:

ec_solution = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = total_enzyme_capacity,
    settings = [set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
    #unconstrain_reactions = ["EX_glc__D_e"],
    optimizer = Tulip.Optimizer,
)

#src these values should be unique (glucose transporter is the only way to get carbon into the system)
@test isapprox(ec_solution.objective, 1.671357282901553, atol = TEST_TOLERANCE) #src
@test isapprox(ec_solution.total_proteome_bound, 0.1, atol = TEST_TOLERANCE) #src
@test isapprox(ec_solution.fluxes.EX_glc__D_e, -49.92966287110028, atol = 0.1) #src
@test isapprox(ec_solution.enzymes.b2417, 0.00011859224858442563, atol = 1e-7) #src

### Building a model incrementally

# Sometimes it is necessary to build a more complicated model, perhaps using a
# novel type of constraint. For this, it is useful to build the enzyme
# constrained model incrementally, using the ConstraintTree building blocks.

import ConstraintTrees as C

# create basic flux model
m = flux_balance_constraints(model)

# create enzyme variables
m += :enzymes^gene_product_variables(model)

# constrain some fluxes...
m.fluxes.EX_glc__D_e.bound = C.Between(-1000.0, 0.0) # undo glucose important bound from original model

# ...And enzymes manually
m.enzymes.b2417.bound = C.Between(0.0, 0.1) # for fun, change the bounds of the protein b2417

# attach the enzyme mass balances
m = with_enzyme_constraints(
    m,
    reaction_isozymes;
    fluxes = m.fluxes, # mount enzyme constraints to these fluxes
    enzymes = m.enzymes, # enzyme variables
)

# add capacity limitation
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        gene_molar_masses,
        A.genes(model),
        total_enzyme_capacity,
    )

# solve the model
ec_solution = optimized_constraints(
    m;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    settings = [set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
