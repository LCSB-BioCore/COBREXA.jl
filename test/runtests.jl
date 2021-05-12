using COBREXA, Test

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
using MCMCChains

# tolerance for comparing analysis results (should be a bit bigger than the
# error tolerance in computations)
TEST_TOLERANCE = 10 * COBREXA._constants.tolerance

function run_test_file(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    @info "$(fn) done in $(round(t; digits = 2))s"
end

function run_test_dir(dir, comment = "Directory $dir/")
    @testset "$comment" begin
        run_test_file.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))
    end
end

function check_data_file_hash(path, expected_checksum)
    actual_checksum = bytes2hex(sha256(open(path)))
    if actual_checksum != expected_checksum
        @error "The downloaded data file `$path' seems to be different from the expected one. Tests will likely fail." actual_checksum expected_checksum
    end
end

function download_data_file(url, path, hash)
    if isfile(path)
        check_data_file_hash(path, hash)
        @info "using cached `$path'"
        return path
    end

    Downloads.download(url, path)
    check_data_file_hash(path, hash)
    return path
end

# set up the workers for Distributed, so that the tests that require more
# workers do not unnecessarily load the stuff multiple times
W = addprocs(2)
@everywhere using COBREXA, Tulip

# load the test models
run_test_file("data", "test_models.jl")

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
end

rmprocs(W)
