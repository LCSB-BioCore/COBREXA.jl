@testset "Small model join" begin
    m1 = load_model(model_paths["e_coli_core.json"])
    m2 = load_model(CoreModel, model_paths["e_coli_core.json"])

    boundary_rxn_ids, boundary_met_ids = boundary_reactions_metabolites(m2)
    exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
    exchange_met_ids = filter(endswith("_e"), boundary_met_ids)

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"]

    community = join_with_exchanges(
        [m1, m2],
        exchange_rxn_ids,
        exchange_met_ids;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(exchange_rxn_ids, reactions(community))
    m2_ex_inds = indexin(exchange_rxn_ids, reactions(m2))
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
    
    boundary_rxn_ids, boundary_met_ids = boundary_reactions_metabolites(m2)
    exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
    exchange_met_ids = filter(endswith("_e"), boundary_met_ids)

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ec_iJO1366_core_53p95M"]

    community = join_with_exchanges(
        [m1, m2],
        exchange_rxn_ids,
        exchange_met_ids;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(exchange_rxn_ids, reactions(community))
    m2_ex_inds = indexin(exchange_rxn_ids, reactions(m2))
    m1_ex_inds = indexin(exchange_rxn_ids, reactions(m1))

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

    boundary_rxn_ids, boundary_met_ids = boundary_reactions_metabolites(m1)
    exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
    exchange_met_ids = filter(endswith("_e"), boundary_met_ids)
    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM"]

    community = join_with_exchanges(
        [m1],
        exchange_rxn_ids,
        exchange_met_ids;
        add_biomass_objective = true,
        biomass_ids = biomass_ids,
    )

    env_ex_inds = indexin(exchange_rxn_ids, reactions(community))
    m1_ex_inds = indexin(exchange_rxn_ids, reactions(m1))
    community.xl[env_ex_inds] .= m1.xl[m1_ex_inds]
    community.xu[env_ex_inds] .= m1.xu[m1_ex_inds]

    m2 = load_model(CoreModel, model_paths["e_coli_core.json"])

    community = add_model_with_exchanges(
        community,
        m2,
        exchange_rxn_ids,
        exchange_met_ids;
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
