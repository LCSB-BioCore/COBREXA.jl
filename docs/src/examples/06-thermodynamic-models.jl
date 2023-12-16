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

flux_solution = Dict(
    "GLCpts" => 1.0,
    "PGI" => 1.0,
    "PFK" => 1.0,
    "FBA" => 1.0,
    "TPI" => 1.0,
    "GAPD" => 1.0,
    "PGK" => -1.0,
    "PGM" => -1.0,
    "ENO" => 1.0,
    "PYK" => 1.0,
    "LDH_D" => -1.0,
)

m = build_max_min_driving_force_model(
    model,
    reaction_standard_gibbs_free_energies;
    flux_solution,
    concentration_lb = 1e-6,
)

m = add_metabolite_ratio_constraints!(
    m,
    Dict("atp" => ("atp_c", "adp_c", 10.0), "nadh" => ("nadh_c", "nad_c", 0.13)),
)

# solve the model
maxmin_solution = optimized_constraints(
    m;
    objective = m.max_min_driving_force.value,
    optimizer = Tulip.Optimizer,
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
