
"""
A named tuple that contains the magic values that are used globally for
whatever purposes.
"""
const _constants = (
    default_reaction_bound = 1e3,
    tolerance = 1e-6,
    sampling_keep_iters = 100,
    sampling_size = 1000,
    possible_rxn_keys = ("rxns", "reactions", "RXNS", "REACTIONS", "Reactions", "Rxns"),
    possible_met_keys = (
        "mets",
        "metabolites",
        "METS",
        "METABOLITES",
        "Metabolites",
        "Mets",
    ),
    possible_gene_keys = ("genes", "GENES", "Genes"),
    possible_lower_bound_keys = ("lbs", "lb", "lowerbounds", "lower_bounds"),
    possible_upper_bound_keys = ("ubs", "ub", "upperbounds", "upper_bounds"),
    possible_stoich_matrix_keys = ("S",),
    possible_balance_keys = ("b",),
    possible_objective_keys = ("c",),
    possible_grr_keys = ("gene_reaction_rules", "grRules", "rules"),
)
