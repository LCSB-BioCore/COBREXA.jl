
# # Flux balance analysis (FBA)

# Here we use [`flux_balance`](@ref) and several related functions to
# find an optimal flux in the *E. coli* "core" model. We will need the model,
# which we can download using [`download_model`](@ref):

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# ## Running a FBA
#
# There are many possibilities on how to arrange the metabolic model into the
# optimization framework and how to actually solve it. The "usual" assumed one
# is captured in the default behavior of function
# [`flux_balance`](@ref):

solution = flux_balance(model, Tulip.Optimizer)

@test isapprox(solution.objective, 0.8739, atol = TEST_TOLERANCE) #src

# The result contains a tree of all optimized values in the model, including
# fluxes, the objective value, and possibly others (given by what the model
# contains).
#
# You can explore the dot notation to explore the solution, extracting e.g. the
# value of the objective:

solution.objective

# ...or the value of the flux through the given reaction (note the solution is
# not unique in FBA):

solution.fluxes.PFK

# ...or make a "table" of all fluxes through all reactions:

collect(solution.fluxes)
