
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
            lattice_args = (samples = 3, ranges = [(-5, 0), (-5, 0), (-5, 5)]),
            workers = W,
        ).values

    @test size(vals) == (3, 3, 3)
    @test count(isnothing, vals) == 15
    @test isapprox(
        filter(!isnothing, vals),
        [
            0.5833914451564178,
            0.5722033617617296,
            0.673404874265479,
            0.5610152783118744,
            0.6606736594947443,
            0.7593405739802911,
            0.7020501075121626,
            0.689318892439569,
            0.787985806919745,
            0.6765876776826879,
            0.7752545924613183,
            0.8739215069575214,
        ],
        atol = TEST_TOLERANCE,
    )
end
