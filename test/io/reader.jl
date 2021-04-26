
@testset "Import MAT" begin
    filepath = joinpath("data", "agora-model.mat")

    download_data_file(
        "https://www.vmh.life/files/reconstructions/AGORA/1.03/reconstructions/mat/Mycoplasma_hominis_ATCC_23114.mat",
        filepath,
        "03362073aa917f0691a0c896948f6e8eebe47f02dcbe0c3f00275fa87396e220",
    )

    cp = load_model(CoreModel, filepath)
    @test cp isa CoreModel
    @test size(cp.S) == (475, 496)
    # @test_throws ErrorException load_model(filepath, "badmodel") # this needs to be fixed

    cp = load_model(CoreModel, joinpath("data", "toyModel1.mat"))
    @test size(cp.S) == (6, 7)

    cp = load_model(CoreModel, joinpath("data", "toyModel2.mat"))
    @test size(cp.S) == (6, 7)

    cp = load_model(CoreModel, joinpath("data", "toyModel3.mat"))
    @test size(cp.S) == (9, 12)

end
