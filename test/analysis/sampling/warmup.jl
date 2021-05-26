@testset "Warm up point generation" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )
    model = load_model(StandardModel, model_path)

    pts = warmup_from_variability(
        100,
        model,
        Tulip.Optimizer;
        modifications = [change_constraint("EX_glc__D_e", -2, 2)],
        workers = W,
    )

    idx = first(indexin(["EX_glc__D_e"], reactions(model)))
    @test size(pts) == (95, 100)
    @test all(pts[idx, :] .>= -2)
    @test all(pts[idx, :] .<= 2)
end
