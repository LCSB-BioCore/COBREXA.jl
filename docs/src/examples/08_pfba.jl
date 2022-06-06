# # Parsimonious flux balance analysis (pFBA)

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use [`flux_balance_analysis`](@ref), [`flux_variability_analysis`](@ref),
# [`parsimonious_flux_balance_analysis`](@ref), and
# [`minimize_metabolic_adjustment_analysis`](@ref), along with the modification functions of
# `COBREXA.jl`, to analyze a toy model of *E. coli*.

# If it is not already present, download the model.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

#md # !!! tip "Tip: use `?` to get quick help about functions"
#md #       When you are unsure about how a function works, write `?
#md #       function_name` to see the function reference documentation.

model = load_model("e_coli_core.xml")

# ## Optimization solvers in `COBREXA`
#
# To actually perform any optimization based analysis we need to load an
# optimizer. Any [`JuMP.jl`-supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here, we will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl)
# to optimize linear programs and
# [`OSQP.jl`](https://osqp.org/docs/get_started/julia.html) to optimize quadratic
# programs.

#md # !!! note "Note: OSQP can be sensitive"
#md #       We recommend reading the docs of `OSQP` before using it, since
#md #       it may give inconsistent results depending on what settings
#md #       you use. Commercial solvers like `Gurobi`, `Mosek`, `CPLEX`, etc.
#md #       require less user engagement.

using Tulip, OSQP, GLPK

# Parsimonious flux balance analysis (here in
# [`parsimonious_flux_balance_analysis`](@ref) finds a unique flux solution
# that minimizes the squared sum of fluxes of the system subject, while
# maintaining the same objective value as the flux balance analysis solution.
# Since we are optimizing a quadratic objective, we also need to switch to a
# quadratic optimizer. In this case, OSQP will work. We demonstrate it on the
# dictionary-returning variant of pFBA,
# [`parsimonious_flux_balance_analysis_dict`](@ref):

dict_soln = parsimonious_flux_balance_analysis_dict(
    model,
    OSQP.Optimizer;
    modifications = [
        silence, # silence the optimizer (OSQP is very verbose by default)
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
        change_optimizer_attribute("polish", true),
    ],
)

# The function also has the expectable second variant that returns a vector of
# solutions, in [`parsimonious_flux_balance_analysis_vec`](@ref). Here, we
# utilize it to show how to use different optimizers for finding the optimum
# and for solving the quadratic problem. That may be preferable if the
# optimizer qualities differ for the differing tasks. pFBA allows you to
# specify `qp_modifications` that are applied after the original optimum is
# found, and before the quadratic part of the problem solving begins.

vec_soln = parsimonious_flux_balance_analysis_vec(
    model,
    Tulip.Optimizer; # start with Tulip
    modifications = [
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
        change_optimizer_attribute("IPM_IterationsLimit", 500), # we may change Tulip-specific attributes here
    ],
    qp_modifications = [
        change_optimizer(OSQP.Optimizer), # now switch to OSQP (Tulip wouldn't be able to finish the computation)
        change_optimizer_attribute("polish", true), # get an accurate solution, see OSQP's documentation
        silence, # and make it quiet.
    ],
)
