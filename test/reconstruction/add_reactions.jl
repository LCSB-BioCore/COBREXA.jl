@testset "@add_reactions! helper" begin
    mod = StandardModel()
    A = Metabolite("A")
    B = Metabolite("B")
    C = Metabolite("C")
    add_metabolites!(mod, A)
    add_metabolites!(mod, B)
    add_metabolites!(mod, C)

    @add_reactions! mod begin
        "v1", nothing ⟷ A
        "v2", nothing ⟷ B, -500
        "v3", nothing ⟷ C, -500, 500
    end

    rxn = mod.reactions["v1"]
    @test rxn.lb == -1000.0
    @test rxn.ub == 1000.0

    rxn = mod.reactions["v2"]
    @test rxn.lb == -500
    @test rxn.ub == 1000.0

    rxn = mod.reactions["v3"]
    @test rxn.lb == -500
    @test rxn.ub == 500
end
