@testset "Summarize" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 200),
        ],
    )
    res = summarize(model, sol; display_unbounded=true, large_flux_bound=20.0)
    @test isnothing(res) # not sure how to test this...
end
