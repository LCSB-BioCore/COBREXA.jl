@testset "Parsimonious flux balance analysis with ObjectModel" begin
    model = test_toyModel()

    d =
        parsimonious_flux_balance_analysis(
            model,
            Tulip.Optimizer;
            modifications = [
                change_constraint("EX_m1(e)", lower_bound = -10.0),
                change_optimizer_attribute("IPM_IterationsLimit", 500),
            ],
            qp_modifications = [change_optimizer(Clarabel.Optimizer), silence],
        ) |> values_dict

    # The used optimizer doesn't really converge to the same answer everytime
    # here, we therefore tolerate a wide range of results.
    @test isapprox(d["biomass1"], 10.0, atol = QP_TEST_TOLERANCE)
end
