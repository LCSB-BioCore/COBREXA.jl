
@testset "Envelopes" begin
    m = load_model(model_paths["e_coli_core.xml"])

    rxns = [1, 2, 3]

    lat = collect.(envelope_lattice(m, rxns; samples = 3))
    @test lat == [[0, 500, 1000], [-1000, 0, 1000], [-1000, 0, 1000]]
    @test lat == collect.(envelope_lattice(m, reactions(m)[rxns]; samples = 3))

    vals =
        objective_envelope(
            m,
            reactions(m)[rxns],
            Tulip.Optimizer;
            lattice = lat,
            workers = W,
        ).values

    @test size(vals) == (3, 3, 3)
    @test count(isnothing, vals) == 15
    @test isapprox(
        filter(!isnothing, vals),
        [
            0.0,
            0.0,
            0.0,
            0.704036947,
            10.053282384,
            10.053282365,
            0.0,
            0.0,
            0.0,
            0.873921506,
            18.664485767,
            20.412058183,
        ],
        atol = TEST_TOLERANCE,
    )
end
