# # Finding balance and variability of constraint-based models

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

# First find a flux distribution that is thermodyamically loopless and incorporates enzyme
# capacity constraints by composing loopless FBA and FBAwMC.

rid_crowding_weight = Dict(# crowding needs a weight for each flux
    rid => 0.004 for rid in reactions(model) if
    !looks_like_biomass_reaction(rid) && !looks_like_exchange_reaction(rid)
)

loopless_crowding_fluxes = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [
        add_crowding_constraints(rid_crowding_weight),
        add_loopless_constraints(),
    ],
)
flux_summary(loopless_crowding_fluxes)