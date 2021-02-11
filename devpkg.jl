using CobraTools
using Measurements
using JLD

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.readmodel(modelpath)

biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
fbasol = CobraTools.fba(model, biomass_rxn)
pfbasol = CobraTools.pfba(model, biomass_rxn)
ad = CobraTools.atom_exchange(pfbasol)

ecoli_kJmol = -71.36 ± 8.55 # formation of biomass kJ/mol

gibbs = CobraTools.mapGibbs(model.rxns) 
JLD.save(joinpath("data", "dgzeros.jld"), "gibbs", gibbs)

CobraTools.blackbox(pfbasol, gibbs, biomass_rxn, ecoli_kJmol)


# # Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(pfbasol)
