@testset "Sampling Tests" begin

    model = load_model(model_paths["e_coli_core.json"])

    cm = MatrixCoupling(model, zeros(1, n_variables(model)), [17.0], [19.0])

    pfk, tala = indexin(["PFK", "TALA"], variables(cm))
    cm.C[:, [pfk, tala]] .= 1.0

    warmup = warmup_from_variability(cm, Tulip.Optimizer; workers = W)

    samples = affine_hit_and_run(
        cm,
        warmup,
        sample_iters = 10 * (1:3),
        workers = W,
        chains = length(W),
    )

    @test size(samples, 1) == size(warmup, 1)
    @test size(samples, 2) == size(warmup, 2) * 3 * length(W)

    lbs, ubs = bounds(model)
    @test all(samples .>= lbs)
    @test all(samples .<= ubs)
    @test all(cm.C * samples .>= cm.cl)
    @test all(cm.C * samples .<= cm.cu)
    @test all(stoichiometry(model) * samples .< TEST_TOLERANCE)
end
