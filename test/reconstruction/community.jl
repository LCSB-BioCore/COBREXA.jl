@testset "Detailed community stoichiometrix matrix check" begin
    m1 = test_toyModel()
    m2 = test_toyModel()
    ex_rxn_mets = Dict("EX_m1(e)" => "m1[e]", "EX_m3(e)" => "m3[e]")

    c1 = join_with_exchanges([m1, m2], ex_rxn_mets; add_biomass_objective = false)

    # test of stoichs are the same
    @test all(c1.S[1:6, 1:7] .== c1.S[7:12, 8:14])
    # test if each models exchange reactions have been added to the environmental exchange properly
    @test sum(c1.S[:, 4]) == 0
    @test sum(c1.S[:, 5]) == 0
    @test sum(c1.S[:, 11]) == 0
    @test sum(c1.S[:, 12]) == 0
    @test sum(c1.S[:, 15]) == -1
    @test sum(c1.S[:, 16]) == -1
    # test if exchange metablites with environment are added properly
    @test c1.S[13, 4] == c1.S[13, 11] == 1
    @test c1.S[14, 5] == c1.S[14, 12] == 1
    # test if environmental exchanges have been added properly
    @test c1.S[13, 15] == c1.S[14, 16] == -1
    # test of bounds set properly
    lb, ub = bounds(c1)
    @test all(lb[1:14] .== -ub[1:14] .== -1000)
    @test all(lb[15:16] .== -ub[15:16] .== 0.0)

    c2 = join_with_exchanges(
        [m1, m2],
        ex_rxn_mets;
        add_biomass_objective = true,
        biomass_ids = ["biomass1", "biomass1"],
    )
    # test if same base stoich matrix
    @test all(c2.S[1:14, 1:16] .== c1.S)
    # test if biomass reaction and metabolites are added correctly
    @test all(c2.S[:, end] .== 0)
    @test c2.S[15, 7] == 1
    @test c2.S[16, 14] == 1

    add_objective!(
        c2,
        ["species_1_biomass1", "species_2_biomass1"];
        objective_weights = [0.1, 0.9],
        objective_column_index = 17,
    )
    @test c2.S[15, end] == -0.1
    @test c2.S[16, end] == -0.9
end

@testset "Small model join" begin
    m1 = load_model(model_paths["e_coli_core.json"])
    m2 = load_model(CoreModel, model_paths["e_coli_core.json"])

    exchange_rxn_mets = Dict(
        ex_rxn => first(keys(reaction_stoichiometry(m2, ex_rxn))) for
        ex_rxn in reactions(m2) if looks_like_exchange_reaction(ex_rxn)
    )

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"]

    community = join_with_exchanges(
        [m1, m2],
        exchange_rxn_mets;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(keys(exchange_rxn_mets), reactions(community))
    m2_ex_inds = indexin(keys(exchange_rxn_mets), reactions(m2))
    community.xl[env_ex_inds] .= m2.xl[m2_ex_inds]
    community.xu[env_ex_inds] .= m2.xu[m2_ex_inds]

    biomass_ids =
        ["species_1_BIOMASS_Ecoli_core_w_GAM", "species_2_BIOMASS_Ecoli_core_w_GAM"]
    add_objective!(
        community,
        biomass_ids;
        objective_column_index = first(
            indexin(["community_biomass"], reactions(community)),
        ),
    )

    d = flux_balance_analysis_dict(community, Tulip.Optimizer)
    @test size(stoichiometry(community)) == (166, 211)
    @test isapprox(d["community_biomass"], 0.41559777495618294, atol = TEST_TOLERANCE)
end

@testset "Heterogenous model join" begin
    m1 = load_model(CoreModel, model_paths["e_coli_core.json"])
    m2 = load_model(CoreModel, model_paths["iJO1366.mat"])

    exchange_rxn_mets = Dict(
        ex_rxn => first(keys(reaction_stoichiometry(m2, ex_rxn))) for
        ex_rxn in reactions(m2) if looks_like_exchange_reaction(ex_rxn)
    )

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ec_iJO1366_core_53p95M"]

    community = join_with_exchanges(
        [m1, m2],
        exchange_rxn_mets;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(keys(exchange_rxn_mets), reactions(community))
    m2_ex_inds = indexin(keys(exchange_rxn_mets), reactions(m2))
    m1_ex_inds = indexin(keys(exchange_rxn_mets), reactions(m1))

    for (env_ex, m2_ex, m1_ex) in zip(env_ex_inds, m2_ex_inds, m1_ex_inds)
        m2lb = isnothing(m2_ex) ? 0.0 : m2.xl[m2_ex]
        m2ub = isnothing(m2_ex) ? 0.0 : m2.xu[m2_ex]

        m1lb = isnothing(m1_ex) ? 0.0 : m1.xl[m1_ex]
        m1ub = isnothing(m1_ex) ? 0.0 : m1.xu[m1_ex]

        community.xl[env_ex] = m1lb + m2lb
        community.xu[env_ex] = m1ub + m2ub
    end

    biomass_metabolite_inds = indexin(
        ["species_1_BIOMASS_Ecoli_core_w_GAM", "species_2_BIOMASS_Ec_iJO1366_core_53p95M"],
        metabolites(community),
    )

    community.S[biomass_metabolite_inds, end] .= -1.0
    community.c[end] = 1.0
    community.xl[end] = 0.0
    community.xu[end] = 1000.0

    d = flux_balance_analysis_dict(
        community,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test size(stoichiometry(community)) == (2203, 3003)
    @test isapprox(d["community_biomass"], 0.8739215069675402, atol = TEST_TOLERANCE)
end

@testset "Community model modifications" begin
    m1 = load_model(CoreModel, model_paths["e_coli_core.json"])

    exchange_rxn_mets = Dict(
        ex_rxn => first(keys(reaction_stoichiometry(m1, ex_rxn))) for
        ex_rxn in reactions(m1) if looks_like_exchange_reaction(ex_rxn)
    )

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM"]

    community = join_with_exchanges(
        [m1],
        exchange_rxn_mets;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(keys(exchange_rxn_mets), reactions(community))
    m1_ex_inds = indexin(keys(exchange_rxn_mets), reactions(m1))
    community.xl[env_ex_inds] .= m1.xl[m1_ex_inds]
    community.xu[env_ex_inds] .= m1.xu[m1_ex_inds]

    m2 = load_model(CoreModel, model_paths["e_coli_core.json"])

    community = add_model_with_exchanges(
        community,
        m2,
        exchange_rxn_mets;
        model_name = "species_2",
        biomass_id = "BIOMASS_Ecoli_core_w_GAM",
    )

    biomass_ids =
        ["species_1_BIOMASS_Ecoli_core_w_GAM", "species_2_BIOMASS_Ecoli_core_w_GAM"]
    add_objective!(
        community,
        biomass_ids;
        objective_column_index = first(
            indexin(["community_biomass"], reactions(community)),
        ),
    )

    d = flux_balance_analysis_dict(community, Tulip.Optimizer)

    @test size(stoichiometry(community)) == (166, 211)
    @test isapprox(d["community_biomass"], 0.41559777495618294, atol = TEST_TOLERANCE)
end
