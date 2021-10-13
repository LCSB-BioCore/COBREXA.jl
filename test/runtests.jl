using COBREXA, Test

using Aqua
using Distributed
using Downloads
using JSON
using JuMP
using LinearAlgebra
using MAT
using OrderedCollections
using OSQP
using SHA
using SparseArrays
using Statistics
using Tulip

# tolerance for comparing analysis results (should be a bit bigger than the
# error tolerance in computations)
TEST_TOLERANCE = 10 * COBREXA._constants.tolerance

print_timing(fn, t) = @info "$(fn) done in $(round(t; digits = 2))s"

# helper functions for running tests en masse
function run_test_file(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    print_timing(fn, t)
end

function run_test_dir(dir, comment = "Directory $dir/")
    @testset "$comment" begin
        run_test_file.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))
    end
end

# set up the workers for Distributed, so that the tests that require more
# workers do not unnecessarily load the stuff multiple times
W = addprocs(2)
t = @elapsed @everywhere using COBREXA, Tulip, JuMP
print_timing("import of packages", t)
t = @elapsed @everywhere begin
    model = Model(Tulip.Optimizer)
    @variable(model, 0 <= x <= 1)
    @objective(model, Max, x)
    optimize!(model)
end
print_timing("JuMP+Tulip code warmup", t)

# make sure there's a directory for temporary data
tmpdir = "tmpfiles"
isdir(tmpdir) || mkdir(tmpdir)
tmpfile(x...) = joinpath(tmpdir, x...)

# load the test models
run_test_file("data_static.jl")
run_test_file("data_downloaded.jl")

# import base files
@testset "COBREXA test suite" begin
    run_test_dir(joinpath("base", "types", "abstract"), "Abstract types")
    run_test_dir(joinpath("base", "types"), "Base model types")
    run_test_dir(joinpath("base", "logging"), "Logging")
    run_test_dir("base", "Base functionality")
    run_test_dir(joinpath("base", "utils"), "Utilities")
    run_test_dir("io", "I/O functions")
    run_test_dir("reconstruction")
    run_test_dir("analysis")
    run_test_dir(joinpath("analysis", "sampling"), "Sampling")
    run_test_file("aqua.jl")
end

rmprocs(W)
