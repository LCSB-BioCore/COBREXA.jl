using Test
using CobraTools

CobraTools.setverbose(false) # quiet

include("testfuncs.jl")

# load model locations
iJO1366_xml = joinpath("..", "models", "iJO1366.xml")
iJO1366_mat = joinpath("..", "models", "iJO1366.mat")
iJO1366_json = joinpath("..", "models", "iJO1366.json") 

yeast_xml = joinpath("..", "models", "yeastGEM.xml")
yeast_mat = joinpath("..", "models", "yeastGEM.mat")

@testset "CobraTools testing" begin
    @testset "Reading & Writing" begin
        # E. coli models
        jsonmodel_ecoli = CobraTools.readmodel(iJO1366_json)
        matlabmodel_ecoli = CobraTools.readmodel(iJO1366_mat)
        sbmlmodel_ecoli = CobraTools.readmodel(iJO1366_xml)

        # Yeast GEM models
        matlabmodel_yeast = CobraTools.readmodel(yeast_mat) # matlab issue
        sbmlmodel_yeast = CobraTools.readmodel(yeast_xml)
        
        @test model_comparison(jsonmodel_ecoli, matlabmodel_ecoli)
        @test_broken model_comparison(sbmlmodel_ecoli, matlabmodel_ecoli) # SBML not implemented yet
        @test_broken model_comparison(matlabmodel_yeast, sbmlmodel_yeast) # yeast GEM matlab format issue

        @test read_write_read(jsonmodel_ecoli, "json")
        @test read_write_read(matlabmodel_ecoli, "mat")
        @test_broken read_write_read(sbmlmodel_ecoli, "xml") # SBML not implemented yet
    end    

    @testset "Construction" begin
        @test_skip false 
    end    

    @testset "Analysis" begin
        @test_skip false 
    end    

    @testset "Sampling" begin
        @test_skip false 
    end

end    

