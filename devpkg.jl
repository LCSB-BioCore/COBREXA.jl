# using Revise

using CobraTools
using MATLAB

jsonpath = joinpath("test/iMM904.json")
jsonmodel = CobraTools.readmodel(jsonpath)

# matpath = joinpath("test/iMM904.mat")
# matpath = joinpath("test", "yeastGEM.mat") # broken


# matmodel = CobraTools.readmodel(matpath)
