# # `MatrixModel` usage

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this tutorial we will introduce `COBREXA`'s `MatrixModel` and
# `MatrixModelWithCoupling`. We will use *E. coli*'s toy model to start with.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

# ## Loading a `MatrixModel`

model = load_model(MatrixModel, "e_coli_core.xml") # we specifically want to load a MatrixModel from the model file

# ## Basic analysis on `MatrixModel`

# As before, for optimization based analysis we need to load an optimizer. Here we
# will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) to optimize the linear
# programs of this tutorial. Refer to the examples of [analysis](05a_fba.md)
# and [analysis modifications](05b_fba_mods.md) for details and explanations.

using Tulip

dict_sol = flux_balance_analysis_dict(
    model,
    Tulip.Optimizer;
    modifications = [
        modify_objective("R_BIOMASS_Ecoli_core_w_GAM"),
        modify_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
        modify_constraint("R_EX_o2_e"; lb = 0, ub = 0),
    ],
)

# ## Structure of `MatrixModel`

# `MatrixModel` is optimized for analysis of models that utilizes the matrix,
# linearly-algebraic "view" of the models. It stores data in a sparse format
# wherever possible.
#
# The structure contains fields that contain the expectable model elements:

fieldnames(MatrixModel)
#
model.S

# Contrary to the usual implementations, the model representation does not
# contain reaction coupling boudns; these can be added to any model by wrapping
# it with [`MatrixCoupling`](@ref). You may also use the prepared
# [`MatrixModelWithCoupling`](@ref) to get a version of [`MatrixModel`](@ref) with this
# coupling.
