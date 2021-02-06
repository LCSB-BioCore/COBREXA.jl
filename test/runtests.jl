using Test

include("testfuncs.jl")

iJO1366_xml = joinpath("test", "iJO1366.xml")
iJO1366_mat = joinpath("test", "iJO1366.mat")
iJO1366_json = joinpath("test", "iJO1366.json") 

iMM904_xml = joinpath("test", "iMM904.xml")
iMM904_mat = joinpath("test", "iMM904.mat")
iMM904_json = joinpath("test", "iMM904.json") 

yeast_xml = joinpath("test", "yeastGEM.xml")
yeast_mat = joinpath("test", "yeastGEM.mat")


    @testset "Test model reading" begin
        @test yeast_xml == yeast_xml 
    end    
    @testset "Test something else" begin
        @test yeast_xml == yeast_xml 
    end    
