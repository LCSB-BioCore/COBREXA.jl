# # Flux balance analysis (FBA)

# We will use [`flux_balance_analysis`](@ref) and several related functions to find the optimal flux in the *E. coli* "core" model.

# If it is not already present, download the model and load the package:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

model = load_model("e_coli_core.xml")

# To perform any optimization based analysis, we need to use a linear problem
# solver (also called an optimizer). Any of the [`JuMP.jl`-supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here, we will demonstrate
# [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) and
# [GLPK](https://www.gnu.org/software/glpk/); other solvers may work just as
# well.

using Tulip

solved_model = flux_balance_analysis(model, Tulip.Optimizer)

# `solved_model` is now an instance of optimized JuMP model. To get the
# variable values out manually, you can use JuMP.value function. Flux variables
# are stored as vector `x`:

using JuMP
value.(solved_model[:x])

# To simplify things, there is a variant of the FBA funtion that does this for
# you automatically:

flux_balance_analysis_vec(model, Tulip.Optimizer)

# Likely, there is another variant that returns the fluxes annotated by
# reaction names, in a dictionary:

flux_balance_analysis_dict(model, Tulip.Optimizer)

# Switching solvers is easy, and may be useful in case you need advanced
# functionality or performance present only in certain solvers. To switch to
# GLPK, you simply load the package and use a different optimizer to run the
# analysis:

using GLPK
flux_balance_analysis_dict(model, GLPK.Optimizer)
