@testset "Construction" begin
    cm = CommunityModel()
    @test isa(cm, CommunityModel)
    cm = CommunityModel(test_toyModel())
    @test isa(cm, CommunityModel)
end

@testset "Basic getters" begin
    cm = CommunityModel(test_toyModel())
    @test reactions(cm) == reactions(test_toyModel())
    @test metabolites(cm) == metabolites(test_toyModel())
    @test stoichiometry(cm) == stoichiometry(test_toyModel())
    cm = CommunityModel(test_LP())
    @test bounds(cm) == bounds(test_LP())
    @test objective(cm) == objective(test_LP())
end

@testset "Actual use case" begin
    m1 = load_model(CoreModel, model_paths["e_coli_core.json"])
    m2 = deepcopy(m1)
    exchange_rxn_mets = Dict(
        ex_rxn => first(keys(reaction_stoichiometry(m2, ex_rxn))) for
        ex_rxn in reactions(m2) if looks_like_exchange_reaction(ex_rxn)
    )
    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"]
    community = join_with_exchanges(
        CoreModel,
        [m1, m2],
        exchange_rxn_mets;
        biomass_ids = biomass_ids,
        model_names = ["m1", "m2"],
    )
    env_ex_inds = indexin(keys(exchange_rxn_mets), reactions(community))
    m2_ex_inds = indexin(keys(exchange_rxn_mets), reactions(m2))
    community.xl[env_ex_inds] .= m2.xl[m2_ex_inds]
    community.xu[env_ex_inds] .= m2.xu[m2_ex_inds]
    biomass_ids =
        Dict("m1_BIOMASS_Ecoli_core_w_GAM" => 1.0, "m2_BIOMASS_Ecoli_core_w_GAM" => 1.0)
    update_community_objective!(community, "community_biomass", biomass_ids)
    cm = CommunityModel(
        community,
        exchange_rxn_mets = exchange_rxn_mets,
        biomass_rxn = "community_biomass",
        model_names = ["m1", "m2"],
    )
    d = flux_balance_analysis_dict(cm, Tulip.Optimizer)
    @test isapprox(d["community_biomass"], 0.41559777495618294, atol = TEST_TOLERANCE)
end
