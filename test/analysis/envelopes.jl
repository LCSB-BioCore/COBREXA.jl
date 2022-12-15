
@testset "Envelopes" begin
    m = load_model(model_paths["e_coli_core.xml"])

    rxns = [1, 2, 3]

    lat = collect.(envelope_lattice(m, rxns; samples = 3))
    @test lat == [[0, 500, 1000], [-1000, 0, 1000], [-1000, 0, 1000]]
    @test lat == collect.(envelope_lattice(m, variables(m)[rxns]; samples = 3))

    vals =
        objective_envelope(
            m,
            variables(m)[rxns],
            Tulip.Optimizer;
            lattice_args = (samples = 3, ranges = [(-5, 5), (-5, 5), (-5, 5)]),
            workers = W,
        ).values

    @test size(vals) == (3, 3, 3)
    @test count(isnothing, vals) == 15
    @test isapprox(
        filter(!isnothing, vals),
        [
            0.0,
            0.0,
            0.3809552834584839,
            0.0,
            0.0,
            0.0,
            0.7619105669322938,
            0.38095528320252253,
            0.0,
            0.0,
            1.0542729305608431,
            0.7619105589042587,
        ],
        atol = TEST_TOLERANCE,
    )
end
