println("loading tests")
@time using Test
@time using Logging
@time using SparseArrays
@time using JuMP
@time using GLPK
@time using Clp
println("loading cobrexa")
@time using COBREXA
println("end cobrexa")
@time using MAT
@time using SHA
@time using Distributed
@time using JuMP
@time using Tulip
@time using OSQP
@time using Statistics
@time using JSON
@time using Measurements
println("end tests")

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

    COBREXA.Downloads.download(url, path)
    check_data_file_hash(path, hash)
    return path
end

include("testing_functions.jl") # load misc. testing functions

# load the test models
run_test_file("data", "test_models.jl")

# import base files
@testset "COBREXA test suite" begin
    run_test_dir("types", "Data structures")
    run_test_dir("base", "Base functionality")
    run_test_dir("io", "I/O functions")
    run_test_dir("reconstruction")
    run_test_dir("analysis")
    run_test_dir("sampling")
end
