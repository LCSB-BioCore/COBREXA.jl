@testset "MOMA" begin
    model = test_toyModel()

    sol = [looks_like_biomass_reaction(rid) ? 0.5 : 0.0 for rid in reactions(model)]

    moma = minimize_metabolic_adjustment_analysis_dict(
        model,
        sol,
        Clarabel.Optimizer;
        modifications = [silence],
    )

    @test isapprox(moma["biomass1"], 0.07692307692307691, atol = QP_TEST_TOLERANCE)
end
