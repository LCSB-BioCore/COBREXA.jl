@testset "Flux variability summary" begin
    model = load_model(model_paths["e_coli_core.json"])

    sol = flux_variability_analysis_dict(
        model,
        Tulip.Optimizer;
        bounds = objective_bounds(0.99),
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 200),
        ],
    )

    fr = flux_variability_summary(sol)
    @test isapprox(fr.biomass_fluxes["BIOMASS_Ecoli_core_w_GAM"][1], 0.8651822872764786; atol=TEST_TOLERANCE)
    @test isapprox(fr.exchange_fluxes["EX_for_e"][2], 1.1440672757160666; atol=TEST_TOLERANCE)
end
