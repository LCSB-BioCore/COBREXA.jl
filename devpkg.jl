using CobraTools

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.readmodel(modelpath)

biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
solobj = CobraTools.fba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(model, biomass_rxn)

# Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# solobj = CobraTools.fba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(model, biomass_rxn)
