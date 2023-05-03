
@testset "Expression-limited models" begin
    orig_model = load_model(model_paths["e_coli_core.json"])

    model =
        orig_model |>
        with_expression_limits(relative_expression = Dict(genes(orig_model) .=> 0.5))

    bs = Dict(reactions(model) .=> bounds(model)[1])

    @test getindex.(Ref(bs), ["PFK", "ENO", "ACALD"]) == [0, -500, -750]

    fluxes = flux_balance_analysis_dict(model, Tulip.Optimizer)

    @test isapprox(fluxes["BIOMASS_Ecoli_core_w_GAM"], 0.21386, atol = TEST_TOLERANCE)
end
