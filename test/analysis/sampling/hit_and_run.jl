@testset "Sampling Tests" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = load_model(StandardModel, model_path)

    using Distributed

    nps = 4 - nprocs()

    # load extra processes, have at least 4 available
    if 1 <= nps <= 3 
        addprocs(nps)
    end
    @everywhere using COBREXA, Tulip

    # Serial test
    chains = hit_and_run(
        model,
        Tulip.Optimizer;
        N = 10_000,
        nchains = 3,
        samplesize = 1000,
        modifications = [
            change_constraint("EX_glc__D_e", -10, -10),
            change_constraint("BIOMASS_Ecoli_core_w_GAM", 0.8, 0.8),
        ],
        workerids = [myid()],
    )

    # Do not test for convergence, that requires too many samples
    @test isapprox(mean(chains[[:PFL]]).nt.mean[1], 0.8, atol = 1.0)

    # parallel tests (still broken, not sure why)
    chainsparallel = hit_and_run(
        model,
        Tulip.Optimizer;
        N = 10_000,
        nchains = 3,
        samplesize = 1000,
        modifications = [
            change_constraint("EX_glc__D_e", -10, -10),
            change_constraint("BIOMASS_Ecoli_core_w_GAM", 0.8, 0.8),
        ],
        workerids = workers(),
    )

    @test isapprox(mean(chainsparallel[[:PFL]]).nt.mean[1], 0.8, atol = 1.0)
end
