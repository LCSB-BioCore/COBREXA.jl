@testset "Parsimonious flux balance analysis with StandardModel" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])
    d = parsimonious_flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint("EX_glc__D_e"; lb = -12, ub = -12),
            change_optimizer_attribute("IPM_IterationsLimit", 500),
        ],
        qp_modifications = [
            change_optimizer(OSQP.Optimizer),
            change_optimizer_attribute("polish", true),
            change_optimizer_attribute("max-iter", 10_000),
            silence,
        ],
    )

    # The used optimizer doesn't really converge to the same answer everytime
    # here, we therefore tolerate a wide range of results.
    @test isapprox(d["PGM"], -17.606459419216442, atol = QP_TEST_TOLERANCE)
end
