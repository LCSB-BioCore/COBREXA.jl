using Revise

using CobraTools

using MATLAB

# jsonpath = joinpath("test/iMM904.json")
# model = CobraTools.readmodel(jsonpath)

matpath = joinpath("test/iMM904.mat")
model = CobraTools.readmodel(jsonpath)
