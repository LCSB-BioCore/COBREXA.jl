@testset "MOMA" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    sol = parsimonious_flux_balance_analysis_dict(
        model,
        OSQP.Optimizer;
        modifications = [silence, change_optimizer_attribute("polish", true)],
    )

    moma = minimize_metabolic_adjustment_analysis_dict(
        model,
        sol,
        OSQP.Optimizer;
        modifications = [
            silence,
            change_optimizer_attribute("polish", true),
            change_constraint("CYTBD"; lb = 0.0, ub = 0.0),
        ],
    )

    @test isapprox(moma["BIOMASS_Ecoli_core_w_GAM"], 0.06214149238730545, atol = 0.05)
end
