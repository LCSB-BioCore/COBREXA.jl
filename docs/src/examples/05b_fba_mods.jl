
# # Extending FBA with modifications

# It is often desirable to add a slight modification to the problem before
# performing the analysis, to see e.g. differences of the model behavior caused
# by the change introduced.
#
# First, let us load everything that will be required:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK, Tulip, JuMP

model = load_model("e_coli_core.xml")

# `COBREXA.jl` supports [many modifications](../concepts/2_modifications.md),
# which include changing objective sense, optimizer attributes, flux
# constraints, optimization objective, reaction and gene knockouts, and others.
# These modifications are applied to the optimization built within the supplied
# optimizer (in this case GLPK) in order as they are specified. User needs to
# manually ensure that the modification ordering is sensible.

# The following example applies multiple different modifications to the *E.
# coli* core model:

fluxes = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [ # modifications are applied in order
        modify_objective("R_BIOMASS_Ecoli_core_w_GAM"), # maximize production
        modify_constraint("R_EX_glc__D_e"; lb = -12, ub = -12), # fix an exchange rate
        knockout(["b0978", "b0734"]), # knock out two genes
        modify_optimizer(Tulip.Optimizer), # ignore the above optimizer and switch to Tulip
        modify_optimizer_attribute("IPM_IterationsLimit", 1000), # customize Tulip
        modify_sense(JuMP.MAX_SENSE), # explicitly tell Tulip to maximize the objective
    ],
)
