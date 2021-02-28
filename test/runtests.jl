using Test
using CobraTools
using JuMP
using Tulip
using OSQP
using Suppressor

include("testfuncs.jl")

# E. coli models
iJO1366_xml = joinpath("..", "models", "iJO1366.xml")
sbmlmodel_ecoli = CobraTools.read_model(iJO1366_xml)

iJO1366_mat = joinpath("..", "models", "iJO1366.mat")
matlabmodel_ecoli = CobraTools.read_model(iJO1366_mat)

iJO1366_json = joinpath("..", "models", "iJO1366.json") 
jsonmodel_ecoli = CobraTools.read_model(iJO1366_json)

core_model = CobraTools.read_model(joinpath("..", "models", "e_coli_core.json")) # to make solving the LPs/QPs faster

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
        @test test_model_manipulations()
    end    

    @testset "Optimization Analysis" begin
        @test fba_test(core_model) 
        @test pfba_test(core_model)  
    end
    
    # @testset "Gibbs Analysis" begin
    #     @test true
    # end

    # @testset "Sampling" begin
    #     @test_skip false 
    # end

end    

