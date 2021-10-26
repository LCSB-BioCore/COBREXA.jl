# # Basic usage of `CoreModel` and `CoreModelCoupled`

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this tutorial we will introduce `COBREXA`'s `CoreModel` and
# `CoreModelCoupled`. We will use *E. coli*'s toy model to start with.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

# ## Loading a `CoreModel`

model = load_model(CoreModel, "e_coli_core.xml") # we specifically want to load a CoreModel from the model file

# ## Basic analysis on `CoreModel`

# As before, for optimization based analysis we need to load an optimizer. Here we
# will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) to optimize the linear
# programs of this tutorial. Refer to the constraint-based analysis basics
# tutorial if you are confused by any functions in this section.

# All the normal analysis functions work on `CoreModel`, due to it also having
# the same generic accessor interface as all the other model types.

using Tulip

dict_sol = flux_balance_analysis_dict(
    model,
    Tulip.Optimizer;
    modifications = [
        change_objective("R_BIOMASS_Ecoli_core_w_GAM"),
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
        change_constraint("R_EX_o2_e"; lb = 0, ub = 0),
    ],
)

# ## Structure of `CoreModel`

# `CoreModel` is a special `COBREXA` type that is optimized for large scale
# analysis of large models. It stores data in a sparse format where possible.

fieldnames(CoreModel)
#
model.S

# ## `CoreModelCoupled` adds coupling constraints to `CoreModel`

# `CoreModelCoupled` extends `CoreModel` by adding coupling constraints.

fieldnames(CoreModelCoupled)

# In short, coupling constraints can be used to ensure that fluxes scale with
# the growth rate (μ) of a model. This reduces the impact of biologically
# infeasible cycles from occurring. Here we will model coupling constraints by
# assuming that they have the form: `-γ ≤ vᵢ/μ  ≤ γ`, where `γ` is the ratio
# between each individual flux (vᵢ) in the model and the growth rate.

gamma = 40 # arbitrary

nr = n_reactions(model) # number of reactions
biomass_index = first(indexin(["R_BIOMASS_Ecoli_core_w_GAM"], reactions(model)))

using LinearAlgebra, SparseArrays

Cf = sparse(1.0I, nr, nr)
Cf[:, biomass_index] .= -gamma
Cb = sparse(1.0I, nr, nr)
Cb[:, biomass_index] .= gamma
C = [Cf; Cb] # coupling constraint matrix

clb = spzeros(2 * nr)
clb[1:nr] .= -1000.0
cub = spzeros(2 * nr)
cub[nr+1:end] .= 1000

cmodel = CoreModelCoupled(model, C, clb, cub)
#
d = flux_balance_analysis_dict(model, Tulip.Optimizer)
d["R_BIOMASS_Ecoli_core_w_GAM"]
#
dc = flux_balance_analysis_dict(cmodel, Tulip.Optimizer)
dc["R_BIOMASS_Ecoli_core_w_GAM"]
