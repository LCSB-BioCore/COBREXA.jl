
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

Load a FBC model representation while guessing the correct model type. Uses
`AbstractFBCModels.load`.

This overload almost always involves a search over types; do not use it in
environments where performance is critical.
"""
function load_model(path::String)
    A.load(path)
end

"""
$(TYPEDSIGNATURES)

Load a FBC model representation. Uses `AbstractFBCModels.load`.
"""
function load_model(model_type::Type{T}, path::String) where {T<:A.AbstractFBCModel}
    A.load(model_type, path)
end

export load_model

"""
$(TYPEDSIGNATURES)

Save a FBC model representation. Uses `AbstractFBCModels.save`.
"""
function save_model(model::T, path::String) where {T<:A.AbstractFBCModel}
    A.save(model, path)
end

export save_model

"""
$(TYPEDSIGNATURES)

Safely download a model with a known hash. All arguments are forwarded to
`AbstractFBCModels.download_data_file` -- see the documentation in the
AbstractFBCModels package for details.
"""
download_model(args...; kwargs...) = A.download_data_file(args...; kwargs...)

export download_model
