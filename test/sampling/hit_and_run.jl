@testset "Sampling Tests" begin
    # these tests are not very good - sampling needs work
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = read_model(model_path, StandardModel)

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    glucose = findfirst(model.reactions, "EX_glc__D_e")
    opt_model = flux_balance_analysis(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective(biomass),
            change_constraint(glucose, -12, -12),
            change_solver_attribute("IPM_IterationsLimit", 500),
        ],
    )
    biomass_index = model[biomass]
    λ = JuMP.value(opt_model[:x][biomass_index])
    change_constraint(biomass, 0.99 * λ, λ)(model, opt_model)

    samples = hit_and_run(100_000, opt_model, keepevery = 10, samplesize = 5000)

    @test isapprox(mean(samples[64, :]), 8.9, atol = 0.1) # only tests if the sampler very approximately converged
end
