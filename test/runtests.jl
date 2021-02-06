using Test
using CobraTools

include("testfuncs.jl")

# load model locations
iJO1366_xml = joinpath("iJO1366.xml")
iJO1366_mat = joinpath("iJO1366.mat")
iJO1366_json = joinpath("iJO1366.json") 

yeast_xml = joinpath("yeastGEM.xml")
yeast_mat = joinpath("yeastGEM.mat")

@testset "Test model reading and writing" begin
    # E. coli models
    jsonmodel_ecoli = CobraTools.readmodel(iJO1366_json)
    matlabmodel_ecoli = CobraTools.readmodel(iJO1366_mat)
    sbmlmodel_ecoli = CobraTools.readmodel(iJO1366_xml)

    # Yeast GEM models
    matlabmodel_yeast = CobraTools.readmodel(yeast_mat)
    sbmlmodel_yeast = CobraTools.readmodel(yeast_xml)
    
    @test model_comparison(jsonmodel_ecoli, matlabmodel_ecoli)
    @test_broken model_comparison(sbmlmodel_ecoli, matlabmodel_ecoli)
    @test_broken model_comparison(matlabmodel_yeast, sbmlmodel_yeast) # yeast GEM matlab format issue

end    

# @testset "Test something else" begin
#     @test yeast_xml == yeast_xml 
# end    
