
@testset "Import MAT model" begin
    cp = load_model(CoreModel, model_paths["mycoplasma-23114.mat"])
    @test cp isa CoreModel
    @test size(cp.S) == (475, 496)
end

@testset "Save MAT model" begin
    testpath = tmpfile("mycoplasma-clone.mat")
    loaded = load_model(CoreModel, model_paths["mycoplasma-23114.mat"])
    save_model(loaded, testpath)
    wrote = load_model(CoreModel, testpath)
    @test wrote isa CoreModel
    @test isequal(wrote, loaded)
end
