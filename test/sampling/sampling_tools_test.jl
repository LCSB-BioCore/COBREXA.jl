@testset "Sampling Tests" begin
    # these tests are not very good - sampling needs work
    model = read_model(
        download_data_file(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
    )

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    glucose = findfirst(model.reactions, "EX_glc__D_e")
    cbm = flux_balance_analysis(model, Tulip.Optimizer; modifications=[modify_objective(biomass), modify_constraint(glucose, -12, -12), modify_solver_attribute("IPM_IterationsLimit", 500)])
    biomass_index = model[biomass]
    λ = JuMP.value(cbm[:x][biomass_index])
    modify_constraint(biomass, 0.99*λ, λ)(model, cbm)

    samples = hit_and_run(
        100_000,
        cbm,
        keepevery = 10,
        samplesize = 5000,
    )

    @test isapprox(mean(samples[64, :]), 8.9, atol = 0.5) # only tests if the sampler approximately converged

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
