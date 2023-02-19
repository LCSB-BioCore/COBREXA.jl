@testset "Flux variability summary" begin
    model = load_model(model_paths["e_coli_core.json"])

    sol =
        flux_variability_analysis_dict(
            model,
            Tulip.Optimizer;
            bounds = objective_bounds(0.90),
            modifications = [change_optimizer_attribute("IPM_IterationsLimit", 2000)],
        ) |> result

    fr = flux_variability_summary(sol)
    @test isapprox(
        fr.biomass_fluxes["BIOMASS_Ecoli_core_w_GAM"][1],
        0.7865293520891825;
        atol = TEST_TOLERANCE,
    )
    @test isapprox(
        fr.exchange_fluxes["EX_for_e"][2],
        11.322324494491848;
        atol = TEST_TOLERANCE,
    )
end
