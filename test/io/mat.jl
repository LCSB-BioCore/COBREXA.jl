
@testset "Import MAT model" begin
    cp = load_model(MatrixModel, model_paths["iJR904.mat"])
    @test cp isa MatrixModel
    @test size(cp.S) == (761, 1075)
end

@testset "Save MAT model" begin
    loaded = load_model(MatrixModel, model_paths["iJR904.mat"])
    testpath = tmpfile("iJR904-clone.mat")
    save_model(loaded, testpath)
    wrote = load_model(MatrixModel, testpath)
    @test wrote isa MatrixModel
    @test isequal(wrote, loaded)
end
