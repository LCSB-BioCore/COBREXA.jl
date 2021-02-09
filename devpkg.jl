using CobraTools
using JLD
using Measurements

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.readmodel(modelpath)

biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
fbasol = CobraTools.fba(model, biomass_rxn)
pfbasol = CobraTools.pfba(model, biomass_rxn)

ad = CobraTools.atom_exchange(pfbasol)

gibbs = CobraTools.mapGibbs(model.rxns)

d = JLD.load(joinpath("data", "bigg_to_kegg.jld"), "bigg_to_kegg")
gibbs_data = JLD.load(joinpath("data", "gibbs_7.0.jld"), "gibbs")

# # Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(pfbasol)
