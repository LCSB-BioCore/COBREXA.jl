using CobraTools

# using SBML
modelpath = joinpath("models", "iMM904.xml") # doesn' work
# modelpath = joinpath("models", "iJO1366.xml") # doesn't work
# modelpath = joinpath("models", "e_coli_core.xml") # doesn't work
# modelpath = joinpath("models", "Ec_core_flux1.xml") # works

# model = readSBML(modelpath) 


# modelpath = joinpath("models", "iAF1260.xml")
# modelpath = joinpath("models", "yeastGEM.xml")


model = CobraTools.read_model(modelpath)
