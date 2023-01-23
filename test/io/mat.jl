
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

@testset "Import yeast-GEM (mat)" begin
    m = load_model(ObjectModel, model_paths["yeast-GEM.mat"])
    @test n_metabolites(m) == 2744
    @test n_reactions(m) == 4063
    @test n_genes(m) == 1160
end
