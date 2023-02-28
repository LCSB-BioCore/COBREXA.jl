@testset "Parsimonious flux balance analysis with ObjectModel" begin
    model = test_toyModel()

    d =
        parsimonious_flux_balance_analysis(
            model,
            Tulip.Optimizer;
            modifications = [
                modify_constraint("EX_m1(e)", lower_bound = -10.0),
                modify_optimizer_attribute("IPM_IterationsLimit", 500),
            ],
            qp_modifications = [modify_optimizer(Clarabel.Optimizer), silence],
        ) |> values_dict

    # The used optimizer doesn't really converge to the same answer everytime
    # here, we therefore tolerate a wide range of results.
    @test isapprox(d["biomass1"], 10.0, atol = QP_TEST_TOLERANCE)

    d2 =
        model |>
        with_changed_bound("biomass1", lower_bound = 10.0) |>
        with_parsimonious_solution(:reaction) |>
        flux_balance_analysis(Clarabel.Optimizer, modifications = [silence]) |>
        values_dict

    @test all(isapprox(d[k], d2[k], atol = QP_TEST_TOLERANCE) for k in keys(d2))

    Q = objective(model |> with_parsimonious_solution(:reaction))
    @test all(Q[i, i] == -1 for i = 1:7)
end
