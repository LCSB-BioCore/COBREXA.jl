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

function runTestFile(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    @info "$(fn) done in $(round(t; digits = 2))s"
end

# load the test models
include(joinpath("data", "testModels.jl"))

# import base files
@testset "COBREXA test suite" begin
    @testset "Base functionality" begin
        runTestFile.("base", ["types.jl", "solver.jl", "utilities.jl"])
    end
    @testset "I/O" begin
        runTestFile.("io", ["reader.jl", "writer.jl", "sbml.jl"])
    end
    @testset "Reconstruction" begin
        runTestFile.("reconstruction", ["coupling.jl", "modeling.jl"])
    end
    @testset "Analysis" begin
        runTestFile.("analysis", ["fba.jl", "fva.jl"])
    end
end
