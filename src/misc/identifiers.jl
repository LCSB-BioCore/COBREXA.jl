"""
    module Identifiers

This module exports interpretation of terms to classify reactions, metabolites,
genes, etc. If an subject has a matching annotation, then it is assumed that it
is part of the associated class of objects. Where possible, this is done using
the terms in module [`SBOTerms`](@ref).

# Exports
$(EXPORTS)
"""
module Identifiers
using ..ModuleTools
@dse

using ..SBOTerms

const EXCHANGE_REACTIONS = [SBOTerms.EXCHANGE_REACTION]

const TRANSPORT_REACTIONS = [
    SBOTerms.TRANSPORT_REACTION,
    SBOTerms.TRANSCELLULAR_MEMBRANE_INFLUX_REACTION,
    SBOTerms.TRANSCELLULAR_MEMBRANE_EFFLUX_REACTION,
    SBOTerms.TRANSLOCATION_REACTION,
    SBOTerms.COTRANSPORT_REACTION,
    SBOTerms.ANTIPORTER_REACTION,
    SBOTerms.SYMPORTER_REACTION,
    SBOTerms.ACTIVE_TRANSPORT,
    SBOTerms.PASSIVE_TRANSPORT,
    SBOTerms.TRANSPORTER,
]

const METABOLIC_REACTIONS = [SBOTerms.BIOCHEMICAL_REACTION]

const BIOMASS_REACTIONS = [SBOTerms.BIOMASS_PRODUCTION]

const ATP_MAINTENANCE_REACTIONS = [SBOTerms.ATP_MAINTENANCE]

const PSEUDO_REACTIONS = [
    SBOTerms.EXCHANGE_REACTION,
    SBOTerms.DEMAND_REACTION,
    SBOTerms.BIOMASS_PRODUCTION,
    SBOTerms.ATP_MAINTENANCE,
    SBOTerms.PSEUDOREACTION,
    SBOTerms.SINK_REACTION,
]

const SPONTANEOUS_REACTIONS = [SBOTerms.SPONTANEOUS_REACTION]

const METABOLITES = [SBOTerms.SIMPLE_CHEMICAL, SBOTerms.METABOLITE]

const GENES = [SBOTerms.GENE]

const REACTIONS = [
    EXCHANGE_REACTIONS
    TRANSPORT_REACTIONS
    METABOLIC_REACTIONS
    BIOMASS_REACTIONS
    ATP_MAINTENANCE_REACTIONS
    PSEUDO_REACTIONS
    SPONTANEOUS_REACTIONS
]

end