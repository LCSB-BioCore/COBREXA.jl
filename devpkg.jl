using CobraTools
using Measurements
using JLD
using JuMP

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.read_model(modelpath)

biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
fbasol = CobraTools.fba(model, biomass_rxn)
pfbasol = CobraTools.pfba(model, biomass_rxn)
ad = CobraTools.atom_exchange(pfbasol)

ecoli_kJmol = -71.36 ± 8.55 # formation of biomass kJ/mol

# gibbs = CobraTools.map_gibbs_rxns(model.rxns, dgtype="prime") # slow 
# CobraTools.saveGibbs(joinpath("data", "dgprime.jld"), gibbs)
gibbs = CobraTools.load_Gibbs(joinpath("data", "dgzeros.jld")) 

CobraTools.map_gibbs_external(pfbasol, gibbs)

CobraTools.map_gibbs_internal(pfbasol, gibbs)




cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)






# # Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(pfbasol)
