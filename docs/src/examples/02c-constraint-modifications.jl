
# # Making adjustments to the constraint system
#
# In the [previous example about model
# adjustments](02b-model-modifications.md), we noted that some constraint
# systems may be to complex to be changed within the limits of the usual FBC
# model view, and we may require a sharper tool to do the changes we need. This
# example shows how to do that by modifying the constraint systems that are
# generated within COBREXA to represent the metabolic model contents.
#
# ## Background: Model-to-optimizer pipeline
#
# ## Background: Constraint trees
#
# ## Changing the model-to-optimizer pipeline
#
# TODO the stuff below:

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# ## Customizing the model

# We can also modify the model. The most explicit way to do this is
# to make a new constraint tree representation of the model.

import ConstraintTrees as C

ctmodel = fbc_model_constraints(model)

fermentation = ctmodel.fluxes.EX_ac_e.value + ctmodel.fluxes.EX_etoh_e.value

forced_mixed_fermentation =
    ctmodel * :fermentation^C.Constraint(fermentation, (10.0, 1000.0)) # new modified model is created

vt = flux_balance_analysis(
    forced_mixed_fermentation,
    Tulip.Optimizer;
    modifications = [silence],
)

@test isapprox(vt.objective, 0.6337, atol = TEST_TOLERANCE) #src

# Models that cannot be solved return `nothing`. In the example below, the
# underlying model is modified.

ctmodel.fluxes.ATPM.bound = (1000.0, 10000.0) # TODO make mutable

vt = flux_balance_analysis(ctmodel, Tulip.Optimizer; modifications = [silence])

@test isnothing(vt) #src

# Models can also be piped into the analysis functions

ctmodel.fluxes.ATPM.bound = (8.39, 10000.0) # revert
vt = ctmodel |> flux_balance_analysis(Tulip.Optimizer; modifications = [silence])

@test isapprox(vt.objective, 0.8739, atol = TEST_TOLERANCE) #src

# Gene knockouts can be done with ease making use of the piping functionality.
# Here oxidative phosphorylation is knocked out.

vt =
    ctmodel |>
    X.knockout!(["b0979", "b0734"], model) |>
    X.flux_balance_analysis(Tulip.Optimizer; modifications = [X.silence])

@test isapprox(vt.objective, 0.21166, atol = TEST_TOLERANCE) #src
