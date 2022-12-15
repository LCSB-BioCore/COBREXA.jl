@testset "MatrixModelWithCoupling generic interface" begin
    model = load_model(MatrixModelWithCoupling, model_paths["e_coli_core.mat"])

    @test reaction_stoichiometry(model, "EX_ac_e") == Dict("ac_e" => -1)
    @test reaction_stoichiometry(model, 44) == Dict("ac_e" => -1)

end

@testset "Conversion from and to ObjectModel" begin
    cm = load_model(MatrixModelWithCoupling, model_paths["e_coli_core.mat"])

    sm = convert(ObjectModel, cm)
    cm2 = convert(MatrixModelWithCoupling, sm)

    @test Set(variables(cm)) == Set(variables(sm))
    @test Set(variables(cm)) == Set(variables(cm2))

    @test reaction_gene_association(sm, variables(sm)[1]) ==
          reaction_gene_association(cm, variables(sm)[1])
end
