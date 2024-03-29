# # Flux balance analysis (FBA)

# We will use [`flux_balance_analysis`](@ref) and several related functions to find the optimal flux in the *E. coli* "core" model.

# If it is not already present, download the model and load the package:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

model = load_model("e_coli_core.xml")

# To perform any optimization-based analysis, we need to use a linear programming
# solver (also called an optimizer). Any of the [`JuMP.jl`-supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here, we will demonstrate
# [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) and
# [GLPK](https://www.gnu.org/software/glpk/); other solvers will likely work just as
# well.

using Tulip

solved_model = flux_balance_analysis(model, Tulip.Optimizer)

# `solved_model` is now an instance of optimized JuMP model. To get the
# variable values out manually, we can use `JuMP.value` function. Flux variables
# are stored as vector `x`:

using JuMP
value.(solved_model[:x])

# To simplify things, there is a variant of the FBA function that does this for
# us automatically:

flux_balance_analysis_vec(model, Tulip.Optimizer)

# Likewise, there is another variant that returns the fluxes annotated by
# reaction names, in a dictionary:

flux_balance_analysis_dict(model, Tulip.Optimizer)

# Switching solvers is easy, and may be useful in case we need advanced
# functionality or performance present only in certain solvers. To switch to
# GLPK, we simply load the package and use a different optimizer to run the
# analysis:

using GLPK
flux_balance_analysis_dict(model, GLPK.Optimizer)

# To get a shortened but useful overview of what was found in the analysis, you
# can use [`flux_summary`](@ref) function:
flux_summary(flux_balance_analysis_dict(model, GLPK.Optimizer))
