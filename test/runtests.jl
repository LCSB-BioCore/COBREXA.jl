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

function runTestDir(dir, comment = "Directory $dir/")
    @testset "$comment" begin
        runTestFile.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))
    end
end

# load the test models
runTestFile("data", "testModels.jl")

# import base files
@testset "COBREXA test suite" begin
    runTestDir("base", "Base functionality")
    runTestDir("io", "I/O functions")
    runTestDir("reconstruction")
    runTestDir("analysis")
end
