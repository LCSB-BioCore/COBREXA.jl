"""
    module SBOTerms

Several selected SBO terms that are recognized by COBREXA. For the full
ontology, see https://github.com/EBI-BioModels/SBO/blob/master/SBO_OBO.obo.

If an SBO term appears here, it *may* be recognized in a function; if an SBO
term does not appear here, then it is *not* used in any COBREXA function.

Mostly, the SBO terms are now used in module `Identifiers`, where they get
grouped semantically to allow other functions to classify reactions (such as
exchanges vs. biomass vs. normal reactions), metabolites, compartments, etc.

# Exports
$(EXPORTS)
"""
module SBOTerms
using ..ModuleTools
@dse

const FLUX_BALANCE_FRAMEWORK = "SBO:0000624"
const RESOURCE_BALANCE_FRAMEWORK = "SBO:0000692"
const CONSTRAINT_BASED_FRAMEWORK = "SBO:0000693"

const PRODUCT = "SBO:0000011"
const CONCENTRATION_OF_PRODUCT = "SBO:0000512"
const SIDE_PRODUCT = "SBO:0000603"

const SUBSTRATE = "SBO:0000015"
const CONCENTRATION_OF_SUBSTRATE = "SBO:0000515"
const SIDE_SUBSTRATE = "SBO:0000604"

const ENZYME = "SBO:0000014"
const TOTAL_CONCENTRATION_OF_ENZYME = "SBO:0000300"

const TRANSCRIPTION = "SBO:0000183"
const TRANSLATION = "SBO:0000184"

const GENE = "SBO:0000243"
const METABOLITE = "SBO:0000299"
const MACROMOLECULE = "SBO:0000245"
const SIMPLE_CHEMICAL = "SBO:0000247"
const RIBONUCLEIC_ACID = "SBO:0000250"
const DEOXYRIBONUCLEIC_ACID = "SBO:0000251"
const TRANSFER_RNA = "SBO:0000313"
const RIBOSOMAL_RNA = "SBO:0000314"
const MESSENGER_RNA = "SBO:0000278"
const TRANSPORTER = "SBO:0000284"
const PROTEIN_COMPLEX = "SBO:0000297"

const MOLECULAR_MASS = "SBO:0000647"
const CATALYTIC_RATE_CONSTANT = "SBO:0000025" # turnover number synonym
const CAPACITY = "SBO:0000661"
const MICHAELIS_CONSTANT = "SBO:0000027"
const MICHAELIS_CONSTANT_FOR_PRODUCT = "SBO:0000323"
const MICHAELIS_CONSTANT_FOR_SUBSTRATE = "SBO:0000322"
const INHIBITORY_CONSTANT = "SBO:0000261"

const STOICHIOMETRIC_COEFFICIENT = "SBO:0000481"
const AND = "SBO:0000173"
const OR = "SBO:0000174"

const PH = "SBO:0000304"
const IONIC_STRENGTH = "SBO:0000623"

const THERMODYNAMIC_TEMPERATURE = "SBO:0000147"
const STANDARD_GIBBS_FREE_ENERGY_OF_REACTION = "SBO:0000583"
const GIBBS_FREE_ENERGY_OF_REACTION = "SBO:0000617"
const STANDARD_GIBBS_FREE_ENERGY_OF_FORMATION = "SBO:0000582"
const TRANSFORMED_STANDARD_GIBBS_FREE_ENERGY_CHANGE_OF_REACTION = "SBO:0000620"
const TRANSFORMED_GIBBS_FREE_ENERGY_CHANGE_OF_REACTION = "SBO:0000622"
const TRANSFORMED_STANDARD_GIBBS_FREE_ENERGY_OF_FORMATION = "SBO:0000621"

const BIOCHEMICAL_OR_TRANSPORT_REACTION = "SBO:0000167"
const BIOCHEMICAL_REACTION = "SBO:0000176"
const TRANSPORT_REACTION = "SBO:0000655"
const TRANSCELLULAR_MEMBRANE_INFLUX_REACTION = "SBO:0000587"
const TRANSCELLULAR_MEMBRANE_EFFLUX_REACTION = "SBO:0000588"
const FLUX_BOUND = "SBO:0000625"
const DEFAULT_FLUX_BOUND = "SBO:0000626"
const EXCHANGE_REACTION = "SBO:0000627"
const DEMAND_REACTION = "SBO:0000628"
const BIOMASS_PRODUCTION = "SBO:0000629"
const ATP_MAINTENANCE = "SBO:0000630"
const PSEUDOREACTION = "SBO:0000631"
const SINK_REACTION = "SBO:0000632"
const SPONTANEOUS_REACTION = "SBO:0000672"
const TRANSLOCATION_REACTION = "SBO:0000185"
const COTRANSPORT_REACTION = "SBO:0000654"
const ANTIPORTER_REACTION = "SBO:0000660"
const SYMPORTER_REACTION = "SBO:0000659"

const ACTIVE_TRANSPORT = "SBO:0000657"
const PASSIVE_TRANSPORT = "SBO:0000658"

const SUBSYSTEM = "SBO:0000633"

end