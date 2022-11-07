using COBREXA.Reconstruction
using COBREXA.Types
using COBREXA.Accessors
using COBREXA.Analysis

using Tulip

m1 = ObjectModel(id = "Model1")
add_metabolites!(
    m1,
    [
        Metabolite("A"),
        Metabolite("B"),
        Metabolite("Ae"),
        Metabolite("Be"),
        Metabolite("X1"),
    ],
)
add_genes!(m1, [Gene("g1"), Gene("g2"), Gene("g3"), Gene("g4")])
add_reactions!(
    m1,
    [
        Reaction("EX_A", Dict("Ae" => -1), :bidirectional),
        Reaction("r1", Dict("Ae" => -1, "A" => 1), :bidirectional),
        Reaction("r2", Dict("A" => -1, "B" => 1, "X1" => 1), :bidirectional),
        Reaction("r3", Dict("B" => -1, "Be" => 1), :forward),
        Reaction("EX_B", Dict("Be" => -1), :forward),
    ],
)

m2 = ObjectModel(id = "Model2")
add_metabolites!(
    m2,
    [
        Metabolite("Ae"),
        Metabolite("A"),
        Metabolite("C"),
        Metabolite("Ce"),
        Metabolite("X2"),
    ],
)
add_genes!(m2, [Gene("g1"), Gene("g2"), Gene("g3"), Gene("g4")])
add_reactions!(
    m2,
    [
        Reaction("r3", Dict("C" => -1, "Ce" => 1), :forward),
        Reaction("EX_C", Dict("Ce" => -1), :forward),
        Reaction("EX_A", Dict("Ae" => -1), :bidirectional),
        Reaction("r1", Dict("Ae" => -1, "A" => 1), :bidirectional),
        Reaction("r2", Dict("A" => -1, "C" => 1, "X2" => 1), :bidirectional),
    ],
)

cm1 = CommunityMember(
    id = "m1",
    abundance = 0.2,
    model = m1,
    exchange_reaction_ids = ["EX_A", "EX_B"],
    biomass_metabolite_id = "X1",
)
cm2 = CommunityMember(
    id = "m2",
    abundance = 0.8,
    model = m2,
    exchange_reaction_ids = ["EX_A", "EX_C"],
    biomass_metabolite_id = "X2",
)


cm = CommunityModel(members = [cm1, cm2], env_met_flux_bounds = Dict("Ae" => (-10, 10)))

reactions(cm)
metabolites(cm)
genes(cm)

n_reactions(cm)
n_metabolites(cm)
n_genes(cm)

stoichiometry(cm)
bounds(cm)
objective(cm)



d = flux_balance_analysis_dict(cm, Tulip.Optimizer)
