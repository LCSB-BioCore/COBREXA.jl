using Test
using CobraTools

include("testfuncs.jl")

# E. coli models
iJO1366_xml = joinpath("..", "models", "iJO1366.xml")
sbmlmodel_ecoli = CobraTools.read_model(iJO1366_xml)

iJO1366_mat = joinpath("..", "models", "iJO1366.mat")
matlabmodel_ecoli = CobraTools.read_model(iJO1366_mat)

iJO1366_json = joinpath("..", "models", "iJO1366.json") 
jsonmodel_ecoli = CobraTools.read_model(iJO1366_json)

@testset "CobraTools Tests" begin

    @testset "Base" begin
        @test test_gene()
        @test test_metabolite()
        @test test_reaction()
        @test test_model()
    end    

    @testset "IO" begin
        @test length(jsonmodel_ecoli.reactions) == 2583
        @test length(matlabmodel_ecoli.reactions) == 2583
        @test_broken length(sbmlmodel_ecoli.reactions) == 2583
        
        @test model_comparison_test(jsonmodel_ecoli, matlabmodel_ecoli)
        @test_broken model_comparison_test(jsonmodel_ecoli, sbmlmodel_ecoli)

        @test read_write_read_test(jsonmodel_ecoli, "json")
        @test read_write_read_test(matlabmodel_ecoli, "mat")
        @test_broken read_write_read_test(sbmlmodel_ecoli, "xml")    
    end    

    @testset "Construction" begin
        @test rxn_construction_test(jsonmodel_ecoli)
        @test (CobraTools.ismassbalanced(findfirst(jsonmodel_ecoli.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M"))[1] == false) && (CobraTools.ismassbalanced(findfirst(jsonmodel_ecoli.rxns, "APCS"))[1])
    end    

    # @testset "Basic Analysis" begin
    #     @test fba_test(jsonmodel_ecoli) 
    #     @test pfba_test(jsonmodel_ecoli)
    #     @test atom_test(jsonmodel_ecoli)
         
    # end
    
    # @testset "Gibbs Analysis" begin
    #     @test true
    # end

    # @testset "Sampling" begin
    #     @test_skip false 
    # end

end    

