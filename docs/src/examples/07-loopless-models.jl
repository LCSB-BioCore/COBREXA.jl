
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

# Additionally to COBREXA and the model format package, we will need a solver --
# let's use GLPK here because we will need to solve mixed interger linear
# programs (MILPs):

import JSONFBCModels
import GLPK
import AbstractFBCModels as A

model = load_model("e_coli_core.json")

# ## Setting up the loopless model

# ## Running a ll-FBA

sol = loopless_flux_balance_analysis(
    model;
    optimizer = GLPK.Optimizer,
    max_flux_bound = 1000.0, # needs to be an order of magnitude bigger, big M method heuristic
    strict_inequality_tolerance = 1.0, # heuristic from paper
    modifications = [set_optimizer_attribute("tol_int", 1e-9)]
)

sol.pseudo_gibbs_free_energy_reaction

@test isapprox(mmdf_solution.max_min_driving_force, 0.8739215069684292, atol = TEST_TOLERANCE) #src

