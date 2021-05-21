@testset "Warm up point generation" begin
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
    ws, lbs, ubs = warmup(
        model,
        Tulip.Optimizer;
        modifications = [change_constraint("EX_glc__D_e", -4, 4)],
        warmup_points = collect(1:n_reactions(model)),
        workerids = [myid()],
    )

    ind = first(indexin(["EX_glc__D_e"], reactions(model)))
    @test size(ws) == (95, 2)
    @test size(ws[1,1]) == (95,)
    @test lbs[ind] ≈ -4
    @test ubs[ind] ≈ 4
    
    # Parallel test
    wsparallel, lbsparallel, ubsparallel = warmup(
        model,
        Tulip.Optimizer;
        modifications = [change_constraint("EX_glc__D_e", -4, 4)],
        warmup_points = collect(1:n_reactions(model)),
        workerids = [workers()],
    )

    ind = first(indexin(["EX_glc__D_e"], reactions(model)))
    @test size(wsparallel) == (95, 2)
    @test size(wsparallel[1,1]) == (95,)
    @test lbsparallel[ind] ≈ -4
    @test ubsparallel[ind] ≈ 4
end
