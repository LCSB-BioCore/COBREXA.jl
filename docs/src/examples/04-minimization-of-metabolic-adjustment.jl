
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

reference_fluxes = parsimonious_flux_balance_analysis(model, Clarabel.Optimizer).fluxes

# TODO MOMA from here
