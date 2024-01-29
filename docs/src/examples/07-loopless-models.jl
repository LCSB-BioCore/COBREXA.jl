
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

# # Loopless flux balance analysis (ll-FBA)

# Here we wil add loopless constraints to a flux balance model to ensure that
# the resultant solution is thermodynamically consistent. As before, we will use
# the core *E. coli* model, which we can download using
# [`download_model`](@ref):

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the JSON model format package. We will also need a
# solver that can solve mixed interger linear programs like GLPK.

import JSONFBCModels
import GLPK

model = load_model("e_coli_core.json")

# ## Running a loopless FBA (ll-FBA)

# One can directly use `loopless_flux_balance_analysis` to solve an FBA problem
# based on `model` where loopless constraints are added to all fluxes. This is
# the direct approach.

solution = loopless_flux_balance_analysis(model; optimizer = GLPK.Optimizer)

@test isapprox(solution.objective, 0.8739215069684303, atol = TEST_TOLERANCE) #src
