@testset "Small model join" begin

    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    m1 = load_model(model_path)
    m2 = load_model(CoreModel, model_path)

    boundary_rxn_ids, boundary_met_ids = all_boundaries(m2)
    exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
    exchange_met_ids = filter(endswith("_e"), boundary_met_ids)

    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM","BIOMASS_Ecoli_core_w_GAM"]

    community = COBREXA.join([m1, m2], exchange_rxn_ids, exchange_met_ids; add_biomass_objective=true, biomass_ids=biomass_ids)

    env_ex_inds = indexin(exchange_rxn_ids, reactions(community))
    m2_ex_inds = indexin(exchange_rxn_ids, reactions(m2))
    community.xl[env_ex_inds] .= m2.xl[m2_ex_inds]
    community.xu[env_ex_inds] .= m2.xu[m2_ex_inds]

    biomass_metabolite_inds = indexin(["species_1_BIOMASS_Ecoli_core_w_GAM", "species_2_BIOMASS_Ecoli_core_w_GAM"], metabolites(community))

    community.S[biomass_metabolite_inds, end] .= -1.0
    community.c[end] = 1.0
    community.xl[end] = 0.0
    community.xu[end] = 1000.0

    d = flux_balance_analysis_dict(community, Tulip.Optimizer)
    @test size(stoichiometry(community)) == (166, 211)
    @test isapprox(d["community_biomass"], 0.41559777495618294, atol = TEST_TOLERANCE,)
end

@testset "Heterogenous model join" begin
    iJO1366_mat = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    )
    core_json = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )


m1 = load_model(core_json)
m2 = load_model(iJO1366_mat)

boundary_rxn_ids, boundary_met_ids = all_boundaries(m2)
exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
exchange_met_ids = filter(endswith("_e"), boundary_met_ids)

biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ec_iML1515_core_75p37M"]

community = COBREXA.join([m1, m2], exchange_rxn_ids, exchange_met_ids; add_biomass_objective=true, biomass_ids=biomass_ids)

exchange_rxn_ids, exchange_met_ids = all_boundaries(community)

exchange_rxn_ids = filter(startswith("EX_"), exchange_rxn_ids)
exchange_met_ids = filter(endswith("_e"), exchange_met_ids)

env_ex_inds = indexin(exchange_rxn_ids, reactions(community))
m2_ex_inds = indexin(exchange_rxn_ids, reactions(m2))
community.xl[env_ex_inds] .= m2.xl[m2_ex_inds]
community.xu[env_ex_inds] .= m2.xu[m2_ex_inds]

biomass_metabolite_inds = indexin(["Core1_BIOMASS_Ecoli_core_w_GAM", "Core2_BIOMASS_Ecoli_core_w_GAM"], metabolites(community))
biomass_metabolite_inds = indexin(["species_1_BIOMASS_Ecoli_core_w_GAM", "species_2_BIOMASS_Ec_iML1515_core_75p37M"], metabolites(community))
community.S[biomass_metabolite_inds, end] .= -1.0
community.c[end] = 1.0
community.xl[end] = 0.0
community.xu[end] = 1000.0

d = flux_balance_analysis_dict(community, Tulip.Optimizer)#; modifications=[change_optimizer_attribute("OutputFlag", 0)])
d["community_biomass"]

end
