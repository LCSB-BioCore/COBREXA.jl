@testset "Parsimonious flux balance analysis with StandardModel" begin
    model = test_toyModel()

    d = parsimonious_flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint("EX_m1(e)", lb = -10.0),
            change_optimizer_attribute("IPM_IterationsLimit", 500),
        ],
        qp_modifications = [
            change_optimizer(OSQP.Optimizer),
            change_optimizer_attribute("polish", true),
            silence,
        ],
    )

    # The used optimizer doesn't really converge to the same answer everytime
    # here, we therefore tolerate a wide range of results.
    @test isapprox(d["biomass1"], 10.0, atol = QP_TEST_TOLERANCE)
end
