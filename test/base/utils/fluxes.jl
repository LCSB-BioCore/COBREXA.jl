@testset "flux utilities" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    fluxes = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    consuming, producing = metabolite_fluxes(model, fluxes)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738; atol = TEST_TOLERANCE)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818; atol = TEST_TOLERANCE)
end
