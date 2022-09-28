
@testset "Conversion from and to MATLAB model" begin
    filename = model_paths["iJO1366.mat"]

    mm = load_mat_model(filename)
    sm = convert(StandardModel, mm)
    mm2 = convert(MATModel, sm)

    @test Set(reactions(mm)) == Set(reactions(sm))
    @test Set(reactions(mm)) == Set(reactions(mm2))
end

@testset "MATModel generic interface" begin
    model = load_model(model_paths["e_coli_core.mat"])

    @test reaction_stoichiometry(model, "EX_ac_e") == Dict("ac_e" => -1)
    @test reaction_stoichiometry(model, 44) == Dict("ac_e" => -1)
end
