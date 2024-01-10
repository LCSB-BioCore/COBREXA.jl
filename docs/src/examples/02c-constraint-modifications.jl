
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

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
# TODO clean up the stuff below:

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

ctmodel = fbc_flux_balance_constraints(model)

fermentation = ctmodel.fluxes.EX_ac_e.value + ctmodel.fluxes.EX_etoh_e.value

forced_mixed_fermentation =
    ctmodel * :fermentation^C.Constraint(fermentation, (10.0, 1000.0)) # new modified model is created

vt = optimized_constraints(
    forced_mixed_fermentation,
    objective = forced_mixed_fermentation.objective.value,
    optimizer = Tulip.Optimizer,
    settings = [silence],
)

@test isapprox(vt.objective, 0.6337, atol = TEST_TOLERANCE) #src

# Models that cannot be solved return `nothing`. In the example below, the
# underlying model is modified.

ctmodel.fluxes.ATPM.bound = C.Between(1000.0, 10000.0)

#TODO explicitly show here how false sharing looks like

vt = optimized_constraints(
    ctmodel,
    objective = ctmodel.objective.value,
    optimizer = Tulip.Optimizer,
    settings = [silence],
)

@test isnothing(vt) #src

# Models can also be piped into the analysis functions

ctmodel.fluxes.ATPM.bound = C.Between(8.39, 10000.0) # revert
vt = optimized_constraints(
    ctmodel,
    objective = ctmodel.objective.value,
    optimizer = Tulip.Optimizer,
    settings = [silence],
)

@test isapprox(vt.objective, 0.8739, atol = TEST_TOLERANCE) #src
