
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
    Maybe{X}

Type of optional values.
"""
const Maybe{X} = Union{Nothing,X}

"""
$(TYPEDEF)

Abstract isozyme type that stores all kinetic information required for
constraint-based modeling.

"""
abstract type Isozyme end

"""
$(TYPEDEF)

A simple struct storing information about the isozyme composition, including
subunit stoichiometry and turnover numbers.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SimpleIsozyme <: Isozyme
    gene_product_stoichiometry::Dict{String,Float64}
    kcat_forward::Maybe{Float64} = nothing
    kcat_backward::Maybe{Float64} = nothing
end

export SimpleIsozyme

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`SimpleIsozyme`](@ref) that takes a string gene
reaction rule and converts it into the appropriate format. Assumes the
`gene_product_stoichiometry` for each subunit is 1.
"""
SimpleIsozyme(gids::Vector{String}; kcat_forward::Float64, kcat_backward::Float64) =
    SimpleIsozyme(;
        gene_product_stoichiometry = Dict(gid => 1.0 for gid in gids),
        kcat_forward,
        kcat_backward,
    )

"""
$(TYPEDEF)

Representation of a binary bound, i.e. constrain a variable to only take the
value 0 or 1 exclusively. Requires a mixed integer-capable solver for
optimization.
"""
struct Binary <: C.Bound end

export Binary
