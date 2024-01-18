
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

"""
$(TYPEDSIGNATURES)

Execute a function with arguments given by `args` on `workers`.

This is merely a nice shortcut for `Distributed.pmap` running over a
`Distributed.CachingPool` of the given workers.
"""
screen(f, args...; workers = D.workers()) = D.pmap(f, D.CachingPool(workers), args...)

"""
$(TYPEDSIGNATURES)

Execute a function arguments from arrays `args` on `workers`, with a pre-cached
JuMP optimization model created from `constraints`, `objective` and `optimizer`
using [`optimization_model`](@ref). `settings` are applied to the optimization
model before first execution of `f`.

Since the model is cached and never re-created, this may be faster than just
plain [`screen`](@ref) in many use cases.

The function `f` is supposed to take `length(args)+1` arguments, the first
argument is the JuMP model, and the other arguments are taken from `args` as
with `Distributed.pmap`. While the model may be modified in place, one should
take care to avoid modifications that change results of subsequent invocations
of `f`, as that almost always results in data races and irreproducible
executions. Ideally, all modifications of the model should be either manually
reverted in the invocation of `f`, or the future invocations of `f` must be
able to overwrite them.

`f` may use [`optimized_model`](@ref) to extract results easily w.r.t. some
given `ConstraintTree`.
"""
function screen_optimization_model(
    f,
    constraints::C.ConstraintTree,
    args...;
    objective::Union{Nothing,C.Value} = nothing,
    sense = Maximal,
    optimizer,
    settings = [],
    workers = D.workers(),
)
    worker_cache = worker_local_data(constraints) do c
        om = COBREXA.optimization_model(c; objective, sense, optimizer)
        for s in settings
            s(om)
        end
        om
    end

    D.pmap(
        (as...) -> f(get_worker_local_data(worker_cache), as...),
        D.CachingPool(workers),
        args...,
    )
end
