@testset "Loopless FBA" begin

    model = load_model(model_paths["e_coli_core.json"])

    sol =
        flux_balance_analysis(
            model,
            GLPK.Optimizer;
            modifications = [add_loopless_constraints()],
        ) |> values_dict

    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        0.8739215069684292,
        atol = TEST_TOLERANCE,
    )
end
