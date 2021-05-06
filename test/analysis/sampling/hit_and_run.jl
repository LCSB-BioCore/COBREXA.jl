@testset "Sampling Tests" begin
    # # these tests are not very good - sampling needs work
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = load_model(StandardModel, model_path)
    opt_model = flux_balance_analysis(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective("BIOMASS_Ecoli_core_w_GAM"),
            change_constraint("EX_glc__D_e", -12, -12),
            change_optimizer_attribute("IPM_IterationsLimit", 500),
        ],
    )
    biomass_index = first(indexin(["BIOMASS_Ecoli_core_w_GAM"], reactions(model)))
    λ = JuMP.value(opt_model[:x][biomass_index])
    change_constraint("BIOMASS_Ecoli_core_w_GAM", 0.99 * λ, λ)(model, opt_model)

    samples = hit_and_run(100_000, opt_model, keepevery = 10, samplesize = 5000)

    # # The sampling converges very randomly and extremely approximately, so only
    # # test a rough result
    @test isapprox(mean(samples[64, :]), 8.9, atol = 0.1)
end
