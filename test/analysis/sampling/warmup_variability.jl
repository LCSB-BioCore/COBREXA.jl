@testset "Warm up point generation" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    pts, lbs, ubs = warmup_from_variability(
        model,
        Tulip.Optimizer,
        100;
        modifications = [change_constraint("EX_glc__D_e"; lb = -2, ub = 2)],
        workers = W,
    )

    idx = first(indexin(["EX_glc__D_e"], reactions(model)))
    @test size(pts) == (95, 100)
    @test all(pts[idx, :] .>= -2)
    @test all(pts[idx, :] .<= 2)
    @test lbs[idx] == -2
    @test ubs[idx] == 2
end
