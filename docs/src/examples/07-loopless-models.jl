
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
import AbstractFBCModels as A

model = load_model("e_coli_core.json")

# ## Running a simple loopless FBA (ll-FBA)

# One can directly use `loopless_flux_balance_analysis` to solve an FBA problem
# based on `model` where loopless constraints are added to all fluxes. This is
# the direct approach. 

sol = loopless_flux_balance_analysis(model; optimizer = GLPK.Optimizer)

@test isapprox(
    sol.objective,
    0.8739215069684303,
    atol = TEST_TOLERANCE,
) #src

@test all(
    v * sol.pseudo_gibbs_free_energy_reaction[k] <= 0
    for (k, v) in sol.fluxes if haskey(sol.pseudo_gibbs_free_energy_reaction, k)
) #src

Dict(
    k => (v,  sol.pseudo_gibbs_free_energy_reaction[k], sol.loopless_binary_variables[k])
for (k, v) in sol.fluxes if abs(v) > 0 && haskey(sol.pseudo_gibbs_free_energy_reaction, k)
)

# ## Building your own loopless model

