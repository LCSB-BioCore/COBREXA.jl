using CobraTools

# using Gurobi
# using Tulip
# using Ipopt
# using GLPK

jsonpath = "test/iJO1366.json"
jsonmodel = CobraTools.readmodel(jsonpath)

# cbmodel = CobraTools.initCBM(jsonmodel)

# biomass_rxn = CobraTools.findrxn(jsonmodel, "BIOMASS_Ec_iJO1366_core_53p95M")
# biomass_ind = CobraTools.getindex(jsonmodel, biomass_rxn)

# CobraTools.fba(cbmodel, biomass_ind)
