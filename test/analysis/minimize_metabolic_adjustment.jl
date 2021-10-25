@testset "MOMA" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    sol = parsimonious_flux_balance_analysis_dict(model, OSQP.Optimizer;)

    moma = minimize_metabolic_adjustment_dict(
        model,
        sol,
        OSQP.Optimizer;
        modifications = [change_constraint("CYTBD"; lb = 0.0, ub = 0.0)],
    )

    @test isapprox(
        moma["BIOMASS_Ecoli_core_w_GAM"],
        0.06214149238730545,
        atol = TEST_TOLERANCE,
    )
end
