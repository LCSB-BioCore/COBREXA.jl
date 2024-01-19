
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

# # Thermodynamic models

using COBREXA

# Here we will solve the max min driving force analysis problem using the
# glycolysis pathway of *E. coli*. In essence, the method attempts to find
# metabolite concentrations (NB: not fluxes) that maximize the smallest
# thermodynamic driving force through each reaction. See Noor, et al., "Pathway
#thermodynamics highlights kinetic obstacles in central metabolism.", PLoS
#computational biology, 2014, for more details.

# To do this, we will first need a model that includes glycolysis, which we can
# download if it is not already present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA, and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# ## Thermodynamic data

# We will need ΔᵣG⁰ data for each reaction we want to include in the
# thermodynamic model. To generate this data manually, go to
# https://equilibrator.weizmann.ac.il/. To generate automatically, use the
# eQuilibrator.jl package.

reaction_standard_gibbs_free_energies = Dict{String,Float64}(
    # Units are in kJ/mol
    "ENO" => -3.8108376097261782,
    "FBA" => 23.376920310319235,
    "GAPD" => 0.5307809794271634,
    "GLCpts" => -45.42430981510088,
    "LDH_D" => 20.04059765689044,
    "PFK" => -18.546314942995934,
    "PGI" => 2.6307087407442395,
    "PGK" => 19.57192102020454,
    "PGM" => -4.470553692565886,
    "PYK" => -24.48733600711958,
    "TPI" => 5.621932460512994,
)

# ## Running basic max min driving force analysis

# If a reference flux is not specified, it is assumed that every reaction in the
# model should be included in the thermodynamic model, and that each reaction
# proceeds in the forward direction. This is usually not intended, and can be
# prevented by inputting a reference flux dictionary as shown below. This
# dictionary can be a flux solution, the sign of each flux is used to determine
# if the reaction runs forward or backward.

# ## Using a reference solution

# Frequently it is useful to check the max-min driving force of a specific FBA
# solution. In this case, one is usually only interested in a subset of all the
# reactions in a model. These reactions can be specified as a the
# `reference_flux`, to only compute the MMDF of these reactions, and ignore all
# other reactions.

reference_flux = Dict(
    "ENO" => 1.0,
    "FBA" => 1.0,
    "GAPD" => 1.0,
    "GLCpts" => 1.0,
    "LDH_D" => -1.0,
    "PFK" => 1.0,
    "PGI" => 1.0,
    "PGK" => -1.0,
    "PGM" => -1.0,
    "PYK" => 1.0,
    "TPI" => 1.0,
)

#!!! warning "Only the signs are extracted from the reference solution"
# It is most convenient to pass a flux solution into `reference_flux`, but
# take care to round fluxes near 0 to their correct sign if they should be
# included in the resultant thermodynamic model. Otherwise, remove them from
# reference flux input.

# ## Solving the MMDF problem

mmdf_solution = max_min_driving_force_analysis(
    model,
    reaction_standard_gibbs_free_energies;
    reference_flux,
    concentration_ratios = Dict(
        "atp" => ("atp_c", "adp_c", 10.0),
        "nadh" => ("nadh_c", "nad_c", 0.13),
    ),
    proton_metabolites = ["h_c", "h_e"],
    water_metabolites = ["h2o_c", "h2o_e"],
    concentration_lower_bound = 1e-6, # M
    concentration_upper_bound = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    optimizer = Tulip.Optimizer,
    settings = [set_optimizer_attribute("IPM_IterationsLimit", 1_000)],
)

@test isapprox(mmdf_solution.max_min_driving_force, 5.78353, atol = TEST_TOLERANCE) #src
