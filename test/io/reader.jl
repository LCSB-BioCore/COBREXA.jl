
@testset "Import MAT" begin
    filepath = joinpath("data", "agora-model.mat")

    download_data_file(
        "https://www.vmh.life/files/reconstructions/AGORA/1.03/reconstructions/mat/Mycoplasma_hominis_ATCC_23114.mat",
        filepath,
        "03362073aa917f0691a0c896948f6e8eebe47f02dcbe0c3f00275fa87396e220",
    )

    cp = load_model(filepath, "model")
    @test cp isa LinearModel
    @test size(cp.S) == (475, 496)
    @test_throws ErrorException load_model(filepath, "badmodel")

    cp = load_model(joinpath("data", "toyModel1.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (6, 7)

    cp = load_model(joinpath("data", "toyModel2.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (6, 7)

    cp = load_model(joinpath("data", "toyModel3.mat"), "model")
    @test cp isa LinearModel
    @test size(cp.S) == (9, 12)

end
