
@testset "Import MAT" begin
    filepath = joinpath("data", "agora-model.mat")
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
