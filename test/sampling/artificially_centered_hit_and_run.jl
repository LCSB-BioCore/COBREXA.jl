@testset "Artificially centered hit and run test" begin
    # these tests are not very good - sampling needs work
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = read_model(model_path, StandardModel)

    biomass = model.reactions["BIOMASS_Ecoli_core_w_GAM"]
    glucose = model.reactions["EX_glc__D_e"]
    opt_model = flux_balance_analysis(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective(biomass),
            change_constraint(glucose, -12, -12),
            change_solver_attribute("IPM_IterationsLimit", 500),
        ],
    )
    biomass_index = index_of(biomass, model)
    λ = JuMP.value(opt_model[:x][biomass_index])
    change_constraint(biomass, 0.99 * λ, λ)(model, opt_model)

    # samples = achr(
    #     100_000,
    #     model,
    #     optimizer;
    #     keepevery = 10,
    #     samplesize = 5000,
    #     constraints = cons,
    # )
    # @test isapprox(mean(samples[64, :]), 8.9, atol = 0.5) # only tests if the sampler approximately converged
end
