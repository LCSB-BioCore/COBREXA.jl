
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

Global configuration options for various COBREXA functions, mainly for various
non-interesting function parameters that are too inconvenient to be passed
around manually.

Changing the configuration values at runtime is possible via the global
[`configuration`](@ref) variable.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Configuration

    """
    Prefixes that [`flux_balance_constraints`](@ref) uses for guessing
    which reactions are exchanges.
    """
    exchange_id_prefixes::Vector{String} = ["EX_", "R_EX_"]

    """
    Prefixes that [`flux_balance_constraints`](@ref) uses for guessing
    which reactions are biomass reactions.
    """
    biomass_id_prefixes::Vector{String} = ["BIOMASS_", "R_BIOMASS_"]

    """
    Reaction identifiers that [`flux_balance_constraints`](@ref) considers to
    be ATP maintenance reactions.
    """
    atp_maintenance_ids::Vector{String} = ["ATPM", "R_ATPM"]

    """
    SBO numbers that label exchange reactions for
    [`flux_balance_constraints`](@ref).
    """
    exchange_sbos::Vector{String} = ["SBO:0000627"]

    """
    SBO numbers that label biomass production reactions for
    [`flux_balance_constraints`](@ref).
    """
    biomass_sbos::Vector{String} = ["SBO:0000629"]

    """
    SBO numbers that label ATP maintenance reactions for
    [`flux_balance_constraints`](@ref).
    """
    atp_maintenance_sbos::Vector{String} = ["SBO:0000630"]

    """
    SBO numbers that label metabolite demand reactions for
    [`flux_balance_constraints`](@ref).
    """
    demand_sbos::Vector{String} = ["SBO:0000628"]
end

"""
    const configuration

The configuration object. You can change the contents of configuration to
override the default behavior of some of the functions.

The available options are described by struct [`Configuration`](@ref).
"""
const configuration = Configuration()
