@testset "Sampling Tests" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = load_model(StandardModel, model_path)

    chains = hit_and_run(
        model,
        Tulip.Optimizer;
        N = 5000_000,
        nchains = 3,
        samplesize = 10_000,
        modifications = [change_constraint("EX_glc__D_e",-10, -10), 
                        change_constraint("BIOMASS_Ecoli_core_w_GAM", 0.8, 0.8)]
        )

    # # The sampling converges very slowly, so can't really do an accurate test 
    # with so few samples
    # this test is ugly, there must be a better way to get the mean
    @test isapprox(mean(chains[[:PFL]]).nt.mean[1], 0.8, atol = 0.5)

    # TODO: parallel tests (still broken, not sure why)
end
