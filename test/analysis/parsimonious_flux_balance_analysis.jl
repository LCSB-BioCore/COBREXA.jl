@testset "Parsimonious flux balance analysis with StandardModel" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])
    d = parsimonious_flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint("EX_glc__D_e", -12, -12),
            change_optimizer_attribute("IPM_IterationsLimit", 500),
        ],
        qp_modifications = [change_optimizer(OSQP.Optimizer), silence],
    )
    v = parsimonious_flux_balance_analysis_vec(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint("EX_glc__D_e", -12, -12),
            change_optimizer_attribute("IPM_IterationsLimit", 500),
        ],
        qp_modifications = [change_optimizer(OSQP.Optimizer), silence],
    )

    # The used optimizer doesn't really converge to the same answer everytime
    # here, we therefore tolerate a wide range of results.
    @test isapprox(d["PGM"], -17.568590034769613, atol = 0.5)
    @test isapprox(v[8], -17.568590034769613, atol = 0.5)
end
