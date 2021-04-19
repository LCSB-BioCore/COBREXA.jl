@testset "Export MAT" begin
    filepath = joinpath("data", "toyModel1.mat")
    loaded = read_model(filepath, CoreModel)
    write_model(loaded, "test_model.mat")
    wrote = read_model("test_model.mat", CoreModel)
    @test wrote isa CoreModel
    @test isequal(wrote, loaded)
    rm("test_model.mat")
end
