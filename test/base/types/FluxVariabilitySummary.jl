@testset "Flux variability summary" begin
    model = load_model(model_paths["e_coli_core.json"])

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 200),
        ],
    )

    fr = flux_variability_summary(sol)
end
