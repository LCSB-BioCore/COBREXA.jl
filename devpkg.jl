using CobraTools

modelpath = joinpath("models", "iJO1366.json") 
jsonmodel = CobraTools.readmodel(modelpath)

biomass_rxn = findfirst(jsonmodel.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
solobj = CobraTools.pfba(jsonmodel, biomass_rxn)


# using JuMP
# using Gurobi
# coremodel = CobraTools.CoreModel(jsonmodel)
# cbmodel, v, massbalance, fluxlbs, fluxubs = CobraTools.initCBM(coremodel)

# objective_index = jsonmodel[biomass_rxn]
# @objective(cbmodel, Max, v[objective_index])
# optimize!(cbmodel)

# cbmodel = CobraTools.initCBM(jsonmodel)

# biomass_rxn = CobraTools.findrxn(jsonmodel, "BIOMASS_Ec_iJO1366_core_53p95M")
# biomass_ind = CobraTools.getindex(jsonmodel, biomass_rxn)

# CobraTools.fba(cbmodel, biomass_ind)
