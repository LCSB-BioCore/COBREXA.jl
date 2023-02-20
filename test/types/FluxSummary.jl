@testset "Flux summary" begin
    model = load_model(model_paths["e_coli_core.json"])

    sol =
        flux_balance_analysis(
            model,
            Tulip.Optimizer;
            modifications = [change_optimizer_attribute("IPM_IterationsLimit", 200)],
        ) |> values_dict

    fr = flux_summary(sol; keep_unbounded = true, large_flux_bound = 25)

    @test isapprox(
        fr.biomass_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8739215022690006;
        atol = TEST_TOLERANCE,
    )
    @test isapprox(fr.export_fluxes["EX_co2_e"], 22.80983339307183; atol = TEST_TOLERANCE)
    @test isapprox(fr.import_fluxes["EX_o2_e"], -21.799492758430517; atol = TEST_TOLERANCE)
    @test isapprox(
        fr.unbounded_fluxes["EX_h2o_e"],
        29.175827202663395;
        atol = TEST_TOLERANCE,
    )
end
