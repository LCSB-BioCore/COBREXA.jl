@testset "Export MAT" begin
    filepath = joinpath("data", "toyModel1.mat")
    loaded = load_model(CoreModel, filepath)
    save_model(loaded, "test_model.mat")
    wrote = load_model(CoreModel, "test_model.mat")
    @test wrote isa CoreModel
    @test isequal(wrote, loaded)
end
