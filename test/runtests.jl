using Test
using CobraTools
using JuMP
using Tulip
using OSQP
using Statistics 

const testdir = dirname(@__FILE__)

include("testing_functions.jl") # load some testing functions

tests = ["base/gene_test.jl",
        "base/metabolite_test.jl",
        "base/model_test.jl",
        "base/reaction_test.jl",
        "io/io_tools_test.jl",
        "construction/construction_overloading_test.jl",
        "construction/model_manipulations_test.jl",
        "optimization_analysis/basic_analysis_test.jl",
        "sampling/sampling_tools_test.jl",
        "external/brenda_tests.jl"]

@testset "CobraTools" begin
    for t in tests
        tp = joinpath(testdir, t)
        include(tp)
    end
end