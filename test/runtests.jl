using Test
using CobraTools

CobraTools.setverbose(false) # quiet

include("testfuncs.jl")

# E. coli models
iJO1366_xml = joinpath("..", "models", "iJO1366.xml")
iJO1366_mat = joinpath("..", "models", "iJO1366.mat")
iJO1366_json = joinpath("..", "models", "iJO1366.json") 

jsonmodel_ecoli = CobraTools.readmodel(iJO1366_json)
matlabmodel_ecoli = CobraTools.readmodel(iJO1366_mat)
sbmlmodel_ecoli = CobraTools.readmodel(iJO1366_xml)

# Yeast GEM models
yeast_xml = joinpath("..", "models", "yeastGEM.xml")
yeast_mat = joinpath("..", "models", "yeastGEM.mat")

matlabmodel_yeast = CobraTools.readmodel(yeast_mat) # matlab issue
sbmlmodel_yeast = CobraTools.readmodel(yeast_xml)


@testset "CobraTools Tests" begin
    @testset "Reading & Writing" begin

        # Are the models be read correctly individually?
        @test length(jsonmodel_ecoli.rxns) == 2583
        @test_broken length(matlabmodel_yeast) == 0
        @test_broken length(sbmlmodel_yeast) == 0

        @test model_comparison_test(jsonmodel_ecoli, matlabmodel_ecoli)
        @test_broken model_comparison_test(sbmlmodel_ecoli, matlabmodel_ecoli) # SBML not implemented yet
        @test_broken model_comparison_test(matlabmodel_yeast, sbmlmodel_yeast) # yeast GEM matlab format issue?

        @test read_write_read_test(jsonmodel_ecoli, "json")
        @test read_write_read_test(matlabmodel_ecoli, "mat")
        @test_broken read_write_read_test(sbmlmodel_ecoli, "xml") # SBML not implemented yet
    end    

    @testset "Construction" begin
        @test rxn_construction_test(jsonmodel_ecoli)
        @test (CobraTools.ismassbalanced(findfirst(jsonmodel_ecoli.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M"))[1] == false) && (CobraTools.ismassbalanced(findfirst(jsonmodel_ecoli.rxns, "APCS"))[1])
    end    

    @testset "Basic Analysis" begin
        @test fba_test(jsonmodel_ecoli) 
        @test pfba_test(jsonmodel_ecoli)
        @test atom_test(jsonmodel_ecoli)
         
    end
    
    @testset "Gibbs Analysis" begin
        @test true
    end

    @testset "Sampling" begin
        @test_skip false 
    end

end    

