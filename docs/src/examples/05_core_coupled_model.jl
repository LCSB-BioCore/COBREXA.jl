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

# Most of the generic accessors and modification functions are also implemented
# for both these model types. This includes adding and deleting genes, reactions
# and metabolites. The objective can also be changed, as well as the bounds on
# the reactions. See the StandardModel example section for hints.
