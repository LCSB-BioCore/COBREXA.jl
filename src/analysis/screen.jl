
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

TODO
"""
screen(f, args...; workers = D.workers()) = D.pmap(f, D.CachingPool(workers), args...)

"""
$(TYPEDSIGNATURES)

TODO also point out there's [`optimized_model`](@ref)

"""
function screen_optimization_model(
    f,
    constraints::C.ConstraintTree,
    args...;
    objective::Union{Nothing,C.Value} = nothing,
    optimizer,
    workers = D.workers(),
)
    # TODO maybe settings?
    worker_cache = worker_local_data(constraints) do c
        (c, COBREXA.optimization_model(c; objective, optimizer))
    end

    D.pmap(
        (as...) -> f(get_worker_local_data(worker_cache)..., as...),
        D.CachingPool(workers),
        args...,
    )
end
