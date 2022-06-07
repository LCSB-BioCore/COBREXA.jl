
# # Extending FBA with modifications

# It is often desirable to add a slight modification to the problem before
# performing the analysis, to see e.g. differences of the model behavior caused
# by the change introduced.
#
# First, let us load everything that will be required:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK, Tulip

model = load_model("e_coli_core.xml")

# `COBREXA.jl` supports [many modifications](TODO), which include changing
# objective sense, optimizer attributes, flux constraints, optimization
# objective, reaction and gene knockouts, and others. These modifications are
# applied to the optimization built within the supplied optimizer (in this case
# GLPK) in order as they are specified. User needs to manually ensure that the
# modification ordering is sensible.

# The following example applies multiple different (although partially
# nonsential) modifications to the *E.  Coli* core model: 

fluxes = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [ # modifications are applied in order
        # this changes the objective to maximize the biomass production
        change_objective("R_BIOMASS_Ecoli_core_w_GAM"),

        # this fixes a specific rate of the glucose exchange
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),

        # this knocks out two genes, i.e. constrains their associated reactions to zero.
        knockout(["b0978", "b0734"]), ## the gene IDs are cytochrome oxidase (CYTBD)

        # ignore the optimizer specified above and change it to Tulip
        change_optimizer(Tulip.Optimizer),

        # set a custom attribute of the Tulip optimizer (see Tulip docs for more possibilities)
        change_optimizer_attribute("IPM_IterationsLimit", 110),

        # explicitly tell the optimizer to maximize the new objective
        change_sense(MAX_SENSE),
    ],
)
