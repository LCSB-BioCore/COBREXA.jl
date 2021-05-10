@testset "Construction overloading" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        joinpath("data", "iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    )

    model = load_model(StandardModel, model_path)

    rxn_original = model.reactions["NADH16pp"]
    nadh = model.metabolites["nadh_c"]
    h_c = model.metabolites["h_c"]
    q8 = model.metabolites["q8_c"]
    q8h2 = model.metabolites["q8h2_c"]
    nad = model.metabolites["nad_c"]
    h_p = model.metabolites["h_p"]

    rxn = nadh + 4.0 * h_c + 1.0 * q8 ⟶ 1.0 * q8h2 + 1.0 * nad + 3.0 * h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0

    rxn = 1.0 * nadh + 4.0 * h_c + q8 ← 1.0 * q8h2 + 1.0 * nad + 3.0 * h_p
    @test rxn.lb < 0.0 && rxn.ub == 0.0

    rxn = 1.0 * nadh + 4.0 * h_c + 1.0 * q8 ↔ q8h2 + nad + 3.0 * h_p
    @test rxn.lb < 0.0 && rxn.ub > 0.0

    rxn = 1.0 * nadh → nothing
    @test length(rxn.metabolites) == 1

    rxn = nothing → nadh
    @test length(rxn.metabolites) == 1

    rxn = nothing → 1.0nadh
    @test length(rxn.metabolites) == 1

    rxn = 1.0 * nadh + 4.0 * h_c + 1.0 * q8 ⟶ 1.0 * q8h2 + 1.0 * nad + 3.0 * h_p
    @test prod(values(rxn.metabolites)) == -12
    @test ("q8h2_c" in [x for x in keys(rxn.metabolites)])

    rxn = nadh + 4.0 * h_c + 1.0 * q8 ⟶ 1.0 * q8h2 + 1.0 * nad + 3.0 * h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0

    @test length(h_p + h_p) == 2
    @test length(h_p + h_p + h_p) == 3

    @test length(+([2.0q8, 1.0nad], nadh)) == 3
end
