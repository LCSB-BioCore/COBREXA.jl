
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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

function run_doc_examples()
    for dir in filter(endswith(".jl"), readdir("../docs/src/examples", join = true))
        @testset "docs/$(basename(dir))" begin
            run_test_file(dir)
        end
    end
end

# set up the workers for Distributed, so that the tests that require more
# workers do not unnecessarily load the stuff multiple times
W = addprocs(2)
t = @elapsed @everywhere begin
    using COBREXA
    import Tulip, JuMP
end

# load the test models
run_test_file("data_static.jl")
run_test_file("data_downloaded.jl")

# TODO data_static and data_downloaded need to be interned into the demos.
# Instead let's make a single "doc running directory" that runs all the
# documentation, which doesn't get erased to improve the test caching.

@testset "COBREXA test suite" begin
    run_doc_examples()
    run_test_file("aqua.jl")
end

rmprocs(W)
