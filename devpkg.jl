using CobraTools
using MATLAB
using JSON

# jsonpath = joinpath("test/iMM904.json")
# m = JSON.parsefile(jsonpath)
# jsonmodel = CobraTools.readmodel(jsonpath)

matpath = joinpath("test/iMM904.mat")
# matpath = joinpath("test", "yeastGEM.mat") # broken


matmodel = CobraTools.readmodel(matpath)

# mf = MatFile(matpath)
# model_name = variable_names(mf)[1] # assume model name is the only variable
# modeldict = get_variable(mf, model_name)
# close(mf)