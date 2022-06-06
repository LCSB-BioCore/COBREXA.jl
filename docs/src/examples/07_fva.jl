# # Flux variability analysis

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

using Tulip

# The default FVA in [`flux_variability_analysis`](@ref) returns maximized and
# minimized reaction fluxes in a matrix. Here we use the dictionary variant in
# flux_variability_analysis_dict, to show how to easily access specific fluxes
# from its results.

fva_mins, fva_maxs = flux_variability_analysis_dict(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.99), # the objective function is allowed to vary by ~1% from the FBA optimum
    modifications = [
        change_optimizer_attribute("IPM_IterationsLimit", 500),
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
)

#
fva_maxs["R_EX_ac_e"]["R_EX_ac_e"] # get the maximal acetate exchange flux

# Another option is to display this information using `flux_variability_summary`. This
# pretty printing only works on flux variability analysis results where dictionary keys indicate
# which flux is optimized and the associated value is a flux dictionary.
flux_variability_summary((fva_mins, fva_maxs))

# More sophisticated variants of [`flux_variability_analysis`](@ref) can be used to extract
# specific pieces of information from the solved optimization problems. Here the objective
# value of the minimized flux and the associated biomass growth rate is returned instead
# of every flux.

biomass_idx = first(indexin(["R_BIOMASS_Ecoli_core_w_GAM"], reactions(model))) # index of biomass function
vs = flux_variability_analysis(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.50), # biomass can vary up to 50% less than optimum
    modifications = [
        change_optimizer_attribute("IPM_IterationsLimit", 500),
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
    ret = m ->
        (COBREXA.JuMP.objective_value(m), COBREXA.JuMP.value(m[:x][biomass_idx])), # m is the model and m[:x] extracts the fluxes from the model
)
#
fva_mins = Dict(rxn => flux for (rxn, flux) in zip(reactions(model), vs[:, 1]))
