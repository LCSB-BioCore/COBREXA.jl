using CobraTools
# hello






# modelpath = joinpath("models", "iJO1366.xml")
# modelpath = joinpath("models", "iJO1366.json")
# modelpath = joinpath("models", "iJO1366.mat")
# modelpath = joinpath("models", "yeastGEM.mat")

# model = CobraTools.read_model(modelpath)

# save_model(model, "test.mat")

# model2 = CobraTools.read_model("test.mat")


# using MAT

# my = matread(modelpath)
# m2 = matread("test.mat")


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

