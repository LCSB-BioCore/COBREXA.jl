@testset "Sampling Tests" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = load_model(StandardModel, model_path)

    warmup, lbs, ubs = warmup_from_variability(model, Tulip.Optimizer; workers = W)

    samples = affine_hit_and_run(
        warmup,
        lbs,
        ubs;
        sample_iters = 10 * (1:3),
        workers = W,
        chains = length(W),
    )

    @test size(samples, 1) == size(warmup, 1)
    @test size(samples, 2) == size(warmup, 2) * 3 * length(W)

    @test all(samples .>= lbs)
    @test all(samples .<= ubs)
    @test all(stoichiometry(model) * samples .< TEST_TOLERANCE)
end
