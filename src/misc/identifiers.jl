"""
This module uses annotation identifiers to classify reactions, metabolites,
genes, etc. If an subject has a matching annotation, then it is assumed that it
is part of the associated class of objects.
"""
module Identifiers
using ..SBOTerms

const EXCHANGE_REACTIONS = [SBOTerms.EXCHANGE_REACTION]

const TRANSPORT_REACTIONS = [
    SBOTerms.TRANSPORT_REACTION,
    SBOTerms.TRANSCELLULAR_MEMBRANE_INFLUX_REACTION,
    SBOTerms.TRANSCELLULAR_MEMBRANE_EFFLUX_REACTION,
]

const METABOLIC_REACTIONS = [SBOTerms.BIOCHEMICAL_REACTION]

const BIOMASS_REACTIONS = [SBOTerms.BIOMASS_PRODUCTION]

const ATP_MAINTENANCE_REACTIONS = [SBOTerms.ATP_MAINTENANCE]

const PSEUDOREACTIONS = [
    SBOTerms.EXCHANGE_REACTION,
    SBOTerms.DEMAND_REACTION,
    SBOTerms.BIOMASS_PRODUCTION,
    SBOTerms.ATP_MAINTENANCE,
    SBOTerms.PSEUDOREACTION,
    SBOTerms.SINK_REACTION,
]

const SPONTANEOUS_REACTIONS = [SBOTerms.SPONTANEOUS_REACTION]
end
