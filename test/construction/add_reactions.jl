@testset "@add_reactions! helper" begin
    mod = ObjectModel(id = "testmodel")
    A = Metabolite(id = "A")
    B = Metabolite(id = "B")
    C = Metabolite(id = "C")
    add_metabolites!(mod, [A, B, C])

    @add_reactions! mod begin
        "v1", nothing ↔ A
        "v2", nothing ↔ B, -500
        "v3", nothing ↔ C, -500, 500
    end

    rxn = mod.reactions["v1"]
    @test rxn.lower_bound == -1000.0
    @test rxn.upper_bound == 1000.0

    rxn = mod.reactions["v2"]
    @test rxn.lower_bound == -500
    @test rxn.upper_bound == 1000.0

    rxn = mod.reactions["v3"]
    @test rxn.lower_bound == -500
    @test rxn.upper_bound == 500
end
