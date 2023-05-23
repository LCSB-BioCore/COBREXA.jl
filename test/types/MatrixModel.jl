@testset "MatrixModel generic interface" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.mat"])

    @test reaction_stoichiometry(model, "EX_ac_e") == Dict("ac_e" => -1)
    @test reaction_stoichiometry(model, 44) == Dict("ac_e" => -1)
end

@testset "Conversion from and to ObjectModel" begin
    cm = load_model(MatrixModel, model_paths["e_coli_core.mat"])

    sm = convert(ObjectModel, cm)
    cm2 = convert(MatrixModel, sm)

    @test Set(variables(cm)) == Set(variables(sm))
    @test Set(variables(cm)) == Set(variables(cm2))

    @test sort(sort.(reaction_gene_associations(sm, reactions(sm)[1]))) ==
          sort(sort.(reaction_gene_associations(cm, reactions(sm)[1])))
end