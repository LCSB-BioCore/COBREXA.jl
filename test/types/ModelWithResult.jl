@testset "ModelWithResult" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])
    res = flux_balance_analysis(model, Tulip.Optimizer)
    @test contains(sprint(show, MIME("text/plain"), res), "ModelWithResult")
end
