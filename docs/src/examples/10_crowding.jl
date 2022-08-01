# # FBA with crowding

# Here we will use [`flux_balance_analysis`](@ref) to explore the metabolism of
# the toy *E. coli* model that additionally respects common protein crowding
# constraints. In particular, the model is limited by the amount of protein
# required to run certain reactions. If that data is available, the predictions
# are accordingly more realistic. See [Beg, et al., "Intracellular crowding
# defines the mode and sequence of substrate uptake by Escherichia coli and
# constrains its metabolic activity.", Proceedings of the National Academy of
# Sciences,2007](https://doi.org/10.1073/pnas.0609845104) for more details.
#
# As usual, the same model modification can be transparently used with many
# other analysis functions, including [`flux_variability_analysis`](@ref) and
# [`parsimonious_flux_balance_analysis`](@ref).

# Let's starting with loading the models and packages.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, Tulip

model = load_model("e_coli_core.xml")

# To describe the protein crowding, each of the enzymes that catalyze the
# reactions gets an associated weight per unit of reaction conversion rate. The
# total sum of all weights multiplied by the flux in the model must be lower
# than 1.
#
# The weights are prepared in a dictionary; for simplicity we assume that the
# relative weight of all enzymes is random between 0.002 and 0.005.
# enzymes are of the same size. Reactions that are not present in the
# dictionary (typically exchanges) are ignored.

import Random
Random.seed!(1) # for repeatability of random numbers below

rid_crowding_weight = Dict(
    rid => 0.002 + 0.003 * rand() for rid in reactions(model) if
    !looks_like_biomass_reaction(rid) && !looks_like_exchange_reaction(rid)
)

# With this, the crowding constraints are added with modification
# [`add_crowding_constraints`](@ref):
loopless_crowding_fluxes = flux_balance_analysis_dict(
    model,
    Tulip.Optimizer;
    modifications = [add_crowding_constraints(rid_crowding_weight)],
)
#
flux_summary(loopless_crowding_fluxes)
