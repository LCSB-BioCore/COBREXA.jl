@testset "FBA with crowding constraints" begin
    model = load_model(model_paths["e_coli_core.json"])
    rid_weight = Dict(
        rid => 0.004 for rid in variables(model) if
        !looks_like_biomass_reaction(rid) && !looks_like_exchange_reaction(rid)
    )

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
            add_crowding_constraints(rid_weight),
            change_constraint("EX_glc__D_e"; lower_bound = -6),
        ],
    )

    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        0.491026987015203,
        atol = TEST_TOLERANCE,
    )

    @test isapprox(sol["EX_ac_e"], 0.7084745257320869, atol = TEST_TOLERANCE)
end
