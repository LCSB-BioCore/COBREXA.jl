@testset "FBA with crowding constraints" begin
    model = load_model(model_paths["e_coli_core.json"])
    idxs = find_internal_reactions(model)

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
            add_crowding_constraint(0.004),
            change_constraint("EX_glc__D_e"; lb = -6),
        ],
    )

    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        0.491026987015203,
        atol = TEST_TOLERANCE,
    )

    @test isapprox(sol["EX_ac_e"], 0.7084745257320869, atol = TEST_TOLERANCE)
end
