
@testset "Import MAT" begin
    filepath = joinpath("data", "agora-model.mat")

    downloadDataFile(
        "https://www.vmh.life/files/reconstructions/AGORA/1.03/reconstructions/mat/Mycoplasma_hominis_ATCC_23114.mat",
        filepath,
        "03362073aa917f0691a0c896948f6e8eebe47f02dcbe0c3f00275fa87396e220",
    )

    cp = loadModel(filepath, "model")
    @test cp isa LinearModel
    @test size(cp.S) == (475, 496)
    @test_throws ErrorException loadModel(filepath, "badmodel")

    cp = loadModel(joinpath("data", "toyModel1.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (6, 7)

    cp = loadModel(joinpath("data", "toyModel2.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (6, 7)

    cp = loadModel(joinpath("data", "toyModel3.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (9, 12)

end
