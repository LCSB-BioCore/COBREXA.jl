using CobraTools

# using SBML
# modelpath = joinpath("models", "iMM904.xml") # doesn' work
# modelpath = joinpath("models", "iJO1366.xml") # doesn't work
# modelpath = joinpath("models", "e_coli_core.xml") # doesn't work
# modelpath = joinpath("models", "Ec_core_flux1.xml") # works

# model = readSBML(modelpath) 


# modelpath = joinpath("models", "iAF1260.xml")
# modelpath = joinpath("models", "yeastGEM.xml")

# modelpath = joinpath("models", "e_coli_core.json")
modelpath = joinpath("models", "iJO1366.mat")
# modelpath = joinpath("models", "iJO1366.json")


model = CobraTools.read_model(modelpath)
# using JSON
# m = JSON.parsefile(modelpath)
# atp = Metabolite("atp") 
# adp = Metabolite("adp") 

# anabolism = 10.0 * atp ⟶ 10.0*adp
# anabolism.id = "anabolism"

# mets = [atp]
# rxns = [anabolism]

# model = Model()
# model.id = "Test model"
# add!(model, mets) # missing adp
# add!(model, rxns)

# fix_model!(model) # adp added

