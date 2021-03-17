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
# load some convenience functions for testing
include("testing_functions.jl") 

# import base files
@testset "COBREXA test suite" begin
    for testSet in ["base", "io", "reconstruction", "analysis", "sampling"]
        @testset "$testSet" begin
            runSuite(testSet)
        end
    end
end

