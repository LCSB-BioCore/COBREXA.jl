using CobraTools
using Measurements
# using PyCall

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

ad = CobraTools.atom_exchange(pfbasol)

f = "kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"
CobraTools.getdG0(f, 7.0, "100")

# Q = pyimport("from equilibrator_api import Q_")
# cc = ComponentContribution()

# gibbs = CobraTools.mapGibbs(model.rxns)

# # Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(pfbasol)
