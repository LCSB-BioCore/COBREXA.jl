# # Constraint-based analysis of single cell models

# In this tutorial we will use `COBREXA`'s `flux_balance_analysis`,
# `flux_variability_analysis`, and `parsimonious_flux_balance_analysis`
# functions to analyze a toy model of *E. coli*.

# If it is not already present, load the model.

!isfile("e_coli_core.xml") && download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")
#
using COBREXA
#
model = load_model("e_coli_core.xml")
#
# ## Optimization solvers in `COBREXA`
#
# To actually perform any optimization based analysis we need to load an
# optimizer. Any [`JuMP` supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here we will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) to
# optimize linear programs and
# [OSQP.jl](https://osqp.org/docs/get_started/julia.html) to optimize quadratic
# programs.

using Tulip
using OSQP

lp_optimizer = Tulip.Optimizer
qp_optimizer = OSQP.Optimizer

# ## Flux balance analysis
#
# Most analysis functions come in two flavours that determine how the output
# will be returned. The suffix of the function determines this behaviour, e.g.
# `flux_balance_analysis_vec` returns a vector of fluxes in the order of the
# reactions returned by `reactions` (recall that this is one of the generic
# interface accessors), and `flux_balance_analysis_dict` that returns a
# dictionary mapping reaction ids to fluxes. 

# In both cases there are two required inputs: the `model` and the `optimizer`.

vec_soln = flux_balance_analysis_vec(model, lp_optimizer)
