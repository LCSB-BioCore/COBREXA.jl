# # Enzyme constrained models

using COBREXA

# Here we will construct an enzyme constrained variant of the *E. coli* "core"
# model. We will need the model, which we can download if it is not already present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip

model = load_model("e_coli_core.json")

# Enzyme constrained models require parameters that are usually not used by
# conventional constraint based models. These include reaction specific turnover
# numbers, molar masses of enzymes, and capacity bounds.

import AbstractFBCModels as A

m1 = fbc_model_constraints(model)
m2 = fbc_model_constraints(model)

lbs, ubs = A.bounds(model)
env_ex_rxns = Dict(rid => (lbs[i], ubs[i]) for (i, rid) in enumerate(A.reactions(model)) if startswith(rid, "EX_"))


m = build_community_environment(env_ex_rxns)
m += :bug1^m1
m += :bug2^m2

member_abundances = [(:bug1, 0.2), (:bug2, 0.8)]
m *=
    :environmental_exchange_balances^link_environmental_exchanges(
        m,
        [(:bug1, 0.2), (:bug2, 0.8)],
    )

m *=
    :equal_growth_rate_constraint^equal_growth_rate_constraints(
        [(:bug1, m.bug1.fluxes.:BIOMASS_Ecoli_core_w_GAM.value), (:bug2, m.bug2.fluxes.:BIOMASS_Ecoli_core_w_GAM.value)]
    )

m.bug1.fluxes.EX_glc__D_e.bound = (-1000.0, 1000.0)
m.bug2.fluxes.EX_glc__D_e.bound = (-1000.0, 1000.0)
m.bug1.fluxes.CYTBD.bound = (-10.0, 10.0) # respiration limited

m *= :objective^C.Constraint(m.bug1.fluxes.:BIOMASS_Ecoli_core_w_GAM.value)

using Gurobi

sol = optimized_constraints(m; objective = m.objective.value, optimizer=Gurobi.Optimizer)

Dict(k => v for (k, v) in sol.bug1.fluxes if startswith(string(k), "EX_"))
Dict(k => v for (k, v) in sol.bug2.fluxes if startswith(string(k), "EX_"))

# exchange cytosolic metabolites
mets = [:EX_akg_e, :EX_succ_e, :EX_pyr_e, :EX_acald_e, :EX_fum_e, :EX_mal__L_e]
for met in mets
    m.bug1.fluxes[met].bound = (-1000.0, 1000.0)
    m.bug2.fluxes[met].bound = (-1000.0, 1000.0)
end

sol = optimized_constraints(m; objective = m.objective.value, optimizer=Tulip.Optimizer, modifications=[set_optimizer_attribute("IPM_IterationsLimit", 100000)])

Dict(k => v for (k, v) in sol.bug1.fluxes if startswith(string(k), "EX_"))
Dict(k => v for (k, v) in sol.bug2.fluxes if startswith(string(k), "EX_"))
