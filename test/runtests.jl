using COBREXA, Test
using Aqua

using Clarabel
using Distributed
import AbstractFBCModels as A
using GLPK # for MILPs

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
    run_doc("01-loading-and-saving.jl")
    run_doc("02-flux-balance-analysis.jl")
    run_test_file("aqua.jl")
end

rmprocs(W)
