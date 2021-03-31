@testset "Export MAT" begin
    filepath = joinpath("data", "toyModel1.mat")
    loaded = load_model(filepath, "model")
    write_model("test_model.mat", loaded, "model")
    wrote = load_model("test_model.mat", "model")
    @test wrote isa LinearModel
    @test isequal(wrote, loaded)
    rm("test_model.mat")
end
