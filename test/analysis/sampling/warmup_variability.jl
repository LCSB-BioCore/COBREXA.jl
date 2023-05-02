@testset "Warm up point generation" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    rid = "EX_glc__D_e"
    pts = warmup_from_variability(
        model |> with_changed_bound(rid; lower_bound = -2, upper_bound = 2),
        Tulip.Optimizer,
        100;
        workers = W,
    )

    idx = first(indexin([rid], variable_ids(model)))
    @test size(pts) == (95, 100)
    @test all(pts[idx, :] .>= -2)
    @test all(pts[idx, :] .<= 2)
end
