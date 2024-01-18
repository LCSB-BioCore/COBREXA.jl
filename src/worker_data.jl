
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
$(TYPEDEF)

Helper struct that provides access to local data that are unboxed and cached
directly on distributed workers.

Use with [`get_worker_local_data`](@ref) and `Distributed.CachingPool`.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct worker_local_data
    "The data that is transferred to the remote worker"
    transfer_data::Any
    "The data that is cached on the remote worker"
    local_data::Union{Some,Nothing}
    """
    The function that converts the transferred data to locally-cached data on
    the remote worker
    """
    transform::Function

    """
    $(TYPEDSIGNATURES)

    Conveniently create a data pack to cache on the remote workers. `f`
    receives a single input (the transfer data `x`) and should produce the
    worker-local data.
    """
    worker_local_data(f, x) = new(x, nothing, f)
end

"""
$(TYPEDSIGNATURES)

"Unwrap" the [`worker_local_data`](@ref) on a remote worker to get the
`local_data` out. If required, executes the `transform` function.

Local copies of `transfer_data` are forgotten after the function executes.
"""
function get_worker_local_data(x::worker_local_data)
    if isnothing(x.local_data)
        x.local_data = Some(x.transform(x.transfer_data))
        x.transfer_data = nothing
    end
    some(x.local_data)
end
