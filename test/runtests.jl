using COBREXA, Test
using Aqua

using Clarabel
using Distributed
import AbstractFBCModels as A
using GLPK # for MILPs

# testing constants
const TEST_TOLERANCE = 1e-3
const QP_TEST_TOLERANCE = 1e-2 # for Clarabel

# helper functions for running tests en masse
print_timing(fn, t) = @info "$(fn) done in $(round(t; digits = 2))s"

function run_test_file(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    print_timing(fn, t)
end

function run_doc_ex(path...)
    run_test_file("..", "docs", "src", "examples", path...)
end

# set up the workers for Distributed, so that the tests that require more
# workers do not unnecessarily load the stuff multiple times
W = addprocs(2)
t = @elapsed @everywhere using COBREXA, Tulip, JuMP

# load the test models
run_test_file("data_static.jl")
run_test_file("data_downloaded.jl")

# import base files
@testset "COBREXA test suite" begin
    run_doc_ex("01-loading-and-saving.jl")
    run_doc_ex("02-flux-balance-analysis.jl")
    run_test_file("aqua.jl")
end

rmprocs(W)
