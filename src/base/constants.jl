
"""
A named tuple that contains the magic values that are used globally for
whatever purposes.
"""
const _constants = (
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
    colors = (empty = :dark_gray, payload = :default, key = :cyan),
)
