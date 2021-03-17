@testset "Export MAT" begin
    filepath = joinpath("data", "toyModel1.mat")
    loaded = loadModel(filepath, "model")
    writeModel("testModel.mat", loaded, "model")
    wrote = loadModel("testModel.mat", "model")
    @test wrote isa LinearModel
    @test isequal(wrote, loaded)
    rm("testModel.mat")
end
