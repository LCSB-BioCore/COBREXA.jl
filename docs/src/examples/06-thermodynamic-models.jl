# # Thermodynamic models

using COBREXA

# Here we will solve the max min driving force analysis problem using the *E.
# coli* "core" model. In essence, the method attempts to find metabolite
# concentrations (NB: not fluxes) that maximize the smallest thermodynamic
# driving force through each reaction. The optimization problem solved is:
# ```
# max min -ΔᵣG
# s.t. ΔᵣG = ΔrG⁰ + R T S' ln(C)
#      ΔᵣG ≤ 0
#      ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
# ```
# where `ΔrG` are the Gibbs energies dissipated by the reactions, R is the gas
# constant, T is the temperature, S is the reaction stoichiometry of the model,
# and C is the vector of metabolite concentrations (and their respective lower
# and upper bounds). See Noor, et al., "Pathway thermodynamics highlights
#kinetic obstacles in central metabolism.", PLoS computational biology, 2014.

# To do this, we will need the model, which we can download if it is not already
# present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA, and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

