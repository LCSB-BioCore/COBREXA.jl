using CobraTools
using MATLAB
using JSON

# jsonpath = "test/iMM904.json"
# m = JSON.parsefile(jsonpath)
# jsonmodel = CobraTools.readmodel(jsonpath)
# CobraTools.savemodel(jsonmodel, "test/iMM904_2.json");
# jsonmodel2 = CobraTools.readmodel("test/iMM904_2.json")
# m2 = JSON.parsefile("test/iMM904_2.json")

matpath = joinpath("test/iMM904.mat")
matmodel = CobraTools.readmodel(matpath)
CobraTools.savemodel(matmodel, "test/iMM904_2.mat")
matpath2 = joinpath("test/iMM904_2.mat")
matmodel2 = CobraTools.readmodel(matpath2)


# matpath = joinpath("test", "yeastGEM.mat") # broken

mf = MatFile(matpath2)
model_name = variable_names(mf)[1] # assume model name is the only variable
modeldict = get_variable(mf, model_name)
close(mf)

