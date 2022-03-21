@testset "MatrixModel generic interface" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.mat"])

    @test reaction_stoichiometry(model, "EX_ac_e") == Dict("ac_e" => -1)
    @test reaction_stoichiometry(model, 44) == Dict("ac_e" => -1)
end

@testset "Conversion from and to ObjectModel" begin
    cm = load_model(MatrixModel, model_paths["e_coli_core.mat"])

    sm = convert(ObjectModel, cm)
    cm2 = convert(MatrixModel, sm)

    @test Set(reactions(cm)) == Set(reactions(sm))
    @test Set(reactions(cm)) == Set(reactions(cm2))

    @test reaction_gene_association(sm, reactions(sm)[1]) ==
          reaction_gene_association(cm, reactions(sm)[1])
end