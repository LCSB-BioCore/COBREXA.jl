
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

# # Minimization of metabolic adjustment analysis

# TODO MOMA citation

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

using COBREXA
import AbstractFBCModels.CanonicalModel as CM
import JSONFBCModels
import Clarabel

# TODO this might do the convert immediately as with the old cobrexa...
# probably better have an actual output-type argument tho rather than force the
# guessing.
model = convert(CM.Model, load_model("e_coli_core.json"))

reference_fluxes =
    parsimonious_flux_balance_analysis(
        model,
        Clarabel.Optimizer,
        settings = [silence],
    ).fluxes

# TODO MOMA from here
