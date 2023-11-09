using COBREXA, Test
using Aqua
using Clarabel
using Distributed
using Downloads
using GLPK # for MILPs
using LinearAlgebra
using Serialization
using SHA
using SparseArrays
using Statistics

# tolerance for comparing analysis results (should be a bit bigger than the
# error tolerance in computations)
TEST_TOLERANCE = 10 * COBREXA.Internal.constants.tolerance
QP_TEST_TOLERANCE = 1e-2 # for Clarabel

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

# make sure there's a directory for temporary data
tmpdir = "tmpfiles"
isdir(tmpdir) || mkdir(tmpdir)
tmpfile(x...) = joinpath(tmpdir, x...)

# load the test models
run_test_file("data_static.jl")
run_test_file("data_downloaded.jl")

# import base files
@testset "COBREXA test suite" begin
    run_test_dir("analysis")
    run_test_file("aqua.jl")
end

rmprocs(W)
