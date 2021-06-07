@testset "CoreModelCoupled generic interface" begin
    model = load_model(CoreModelCoupled, model_paths["e_coli_core.mat"])

    @test reaction_equation(model, "EX_ac_e") == Dict("ac_e" => -1)
    @test reaction_equation(model, 44) == Dict("ac_e" => -1)

end
