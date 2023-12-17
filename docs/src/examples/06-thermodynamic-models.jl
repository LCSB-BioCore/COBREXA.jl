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
    "GLCpts" => -45.42430981510088,
    "PGI" => 2.6307087407442395,
    "PFK" => -18.546314942995934,
    "FBA" => 23.376920310319235,
    "TPI" => 5.621932460512994,
    "GAPD" => 0.5307809794271634,
    "PGK" => 19.57192102020454,
    "PGM" => -4.470553692565886,
    "ENO" => -3.8108376097261782,
    "PYK" => -24.48733600711958,
    "LDH_D" => 20.04059765689044,
)

# ## Running basic max min driving force analysis

# If a reference flux is not specified, it is assumed that every reaction in the
# model should be included in the thermodynamic model, and that each reaction
# proceeds in the forward direction. This is usually not intended, and can be
# prevented by inputting a reference flux dictionary as shown below. This
# dictionary can be a flux solution, the sign of each flux is used to determine
# if the reaction runs forward or backward.

#!!! warning "Only the signs are extracted from the reference solution"
# It is most convenient to pass a flux solution into `reference_flux`, but
# take care to round fluxes near 0 to their correct sign if they should be
# included in the resultant thermodynamic model. Otherwise, remove them from
# reference flux input.

reference_flux = Dict(
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

mmdf_solution = max_min_driving_force_analysis(
    model,
    reaction_standard_gibbs_free_energies;
    reference_flux,
    concentration_ratios = Dict(
        "atp" => ("atp_c", "adp_c", 10.0),
        "nadh" => ("nadh_c", "nad_c", 0.13),
    ),
    proton_ids = ["h_c", "h_e"],
    water_ids = ["h2o_c", "h2o_e"],
    concentration_lb = 1e-6, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 1_000)],
    optimizer = Tulip.Optimizer,
)

@test isapprox(mmdf_solution.max_min_driving_force, 5.78353, atol = TEST_TOLERANCE) #src

# ## Building the thermodynamic model yourself

# It is also possible to build the thermodynamic model yourself. This allows you
# to incorporate more complex constraints and gives you more freedom.

m = build_max_min_driving_force_model(
    model,
    reaction_standard_gibbs_free_energies;
    proton_ids = ["h_c", "h_e"],
    water_ids = ["h2o_c", "h2o_e"],
    reference_flux,
    concentration_lb = 1e-6, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
)

m *= :metabolite_ratio_constraints^log_ratio_constraints(
    Dict("atp" => ("atp_c", "adp_c", 10.0), "nadh" => ("nadh_c", "nad_c", 0.13)),
    m.log_metabolite_concentrations,
)

# solve the model
mmdf_solution = optimized_constraints(
    m;
    objective = m.max_min_driving_force.value,
    optimizer = Tulip.Optimizer,
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 1_000)],
)
