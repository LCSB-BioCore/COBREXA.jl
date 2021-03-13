using Test
using Logging
using SparseArrays
using JuMP
using GLPK
using Clp
using COBREXA
using MAT
using SHA
using Distributed
using JuMP

using Tulip
using OSQP
using Statistics
using JSON
using Measurements

"""
    runSuite(baseDir)

Runs the tests located in the `baseDir` and outputs the time it took
"""
function runSuite(baseDir)
    for file in filter(f -> endswith(f, ".jl"), readdir(baseDir))
        t = time()
        include(joinpath(baseDir, file))
        @info "$(file) took $(round(time() - t; digits = 2)) seconds"
    end
end

# load the test models
include(joinpath("data", "testModels.jl"))

# import base files
@testset "COBREXA test suite" begin
    for testSet in ["base", "io", "reconstruction", "analysis"]
        @testset "$testSet" begin
            runSuite(testSet)
        end
    end
end


const testdir = dirname(@__FILE__)

include("testing_functions.jl") # load some testing functions

tests = [
    "base/gene_test.jl",
    "base/metabolite_test.jl",
    "base/model_test.jl",
    "base/reaction_test.jl",
    "io/io_tools_test.jl",
    "construction/construction_overloading_test.jl",
    "construction/model_manipulations_test.jl",
    "optimization_analysis/basic_analysis_test.jl",
    "sampling/sampling_tools_test.jl",
    "external/brenda_tests.jl",
]

@testset "CobraTools" begin
    for t in tests
        tp = joinpath(testdir, t)
        include(tp)
    end
end

