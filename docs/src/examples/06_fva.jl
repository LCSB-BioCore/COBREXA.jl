# # Flux variability analysis (FVA)

# Here we will use [`flux_variability_analysis`](@ref) to analyze the *E. coli*
# core model.

# As usual, it is not already present, download the model and load the required
# packages. We picked the GLPK solver, but others may work as well:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# The FVA implementation in [`flux_variability_analysis`](@ref) returns
# maximized and minimized reaction fluxes in a 2-column matrix.
# The bounds parameter here sets the "gamma" parameter -- the objective
# function is allowed to vary by around 1% from the optimum found by FBA on the
# same model:

flux_variability_analysis(model, GLPK.Optimizer; bounds = objective_bounds(0.99))

# ## Detailed variability analysis with modifications
#
# A dictionary-returning variant in [`flux_variability_analysis_dict`](@ref),
# returns the result in a slightly more structured way. At the same time, we
# can specify additional [modifications](../concepts/2_modifications.md) to be
# applied to the model:

min_fluxes, max_fluxes = flux_variability_analysis_dict(
    model,
    GLPK.Optimizer;
    bounds = objective_bounds(0.99),
    modifications = [
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
)

# The dictionaries can be easily used to explore the whole state of the model
# when certain reactions are maximized or minimized. For example, we can take
# the maximal acetate exchange flux when the acetate exchange is maximized:

max_fluxes["R_EX_ac_e"]["R_EX_ac_e"]

# We can also check that the modifications really had the desired effect on
# oxygen consumption:

max_fluxes["R_EX_ac_e"]["R_O2t"]

# ...and see how much carbon dioxide would produced under at the given
# metabolic extreme:

max_fluxes["R_EX_ac_e"]["R_EX_co2_e"]


# ## Summarizing the flux variability
#
# A convenience function [`flux_variability_summary`](@ref) is able to display
# this information in a nice overview:
flux_variability_summary((min_fluxes, max_fluxes))

# ## Customizing the FVA output
#
# Parameter `ret` of [`flux_variability_analysis`](@ref) can be used to extract
# specific pieces of information from the individual solved (minimized and
# maximized) optimization problems. Here we show how to extract the value of
# biomass "growth" along with the minimized/maximized reaction flux:

# find the index of biomass reaction in all reactions
biomass_idx = first(indexin(["R_BIOMASS_Ecoli_core_w_GAM"], reactions(model)))

vs = flux_variability_analysis(
    model,
    GLPK.Optimizer;
    bounds = objective_bounds(0.50), # objective can vary by up to 50% of the optimum
    modifications = [
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
    ret = optimized_model -> (
        COBREXA.JuMP.objective_value(optimized_model),
        COBREXA.JuMP.value(optimized_model[:x][biomass_idx]),
    ),
)
