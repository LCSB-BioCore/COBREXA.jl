
"""
A named tuple that contains the magic values that are used globally for
whatever purposes.
"""
const _constants = (
    default_stoich_show_size = 50_000,
    default_reaction_bound = 1e3,
    tolerance = 1e-6,
    sampling_keep_iters = 100,
    sampling_size = 1000,
    keynames = (
        rxns = ["reactions", "rxns", "RXNS", "REACTIONS", "Reactions", "Rxns"],
        mets = ["metabolites", "mets", "METS", "METABOLITES", "Metabolites", "Mets"],
        genes = ["genes", "GENES", "Genes"],
        lbs = ["lbs", "lb", "lowerbounds", "lower_bounds"],
        ubs = ["ubs", "ub", "upperbounds", "upper_bounds"],
        stochiometry = ["S"],
        balance = ["b"],
        objective = ["c"],
        grrs = ["gene_reaction_rules", "grRules", "rules"],
        ids = ["id", "description"],
    ),
    gene_annotation_checks = (
        "ncbigene",
        "ncbigi",
        "refseq_locus_tag",
        "refseq_name",
        "refseq_synonym",
        "uniprot",
    ),
    reaction_annotation_checks = (
        "bigg.reaction",
        "biocyc",
        "ec-code",
        "kegg.reaction",
        "metanetx.reaction",
        "rhea",
        "sabiork",
        "seed.reaction",
    ),
    metabolite_annotation_checks = (
        "kegg.compound",
        "bigg.metabolite",
        "chebi",
        "inchi_key",
        "sabiork",
        "hmdb",
        "seed.compound",
        "metanetx.chemical",
        "reactome.compound",
        "biocyc",
    ),
)

const MAX_SENSE = COBREXA.MOI.MAX_SENSE
const MIN_SENSE = COBREXA.MOI.MIN_SENSE
