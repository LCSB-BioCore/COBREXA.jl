@testset "Construction overloading" begin
    model = CobraTools.read_model(joinpath("data", "iJO1366.json") )
    
    rxn_original = findfirst(model.reactions, "NADH16pp")
    nadh = findfirst(model.metabolites, "nadh_c")
    h_c = findfirst(model.metabolites, "h_c")
    q8 = findfirst(model.metabolites, "q8_c")
    q8h2 = findfirst(model.metabolites, "q8h2_c")
    nad = findfirst(model.metabolites, "nad_c")
    h_p = findfirst(model.metabolites, "h_p")
    
    rxn = nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0

    rxn = 1.0*nadh + 4.0*h_c + q8 ← 1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb < 0.0 && rxn.ub == 0.0
  
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ↔ q8h2 + nad + 3.0*h_p
    @test rxn.lb < 0.0 && rxn.ub > 0.0
    
    rxn = 1.0*nadh → ∅
    @test length(rxn.metabolites) == 1

    rxn = ∅ → nadh
    @test length(rxn.metabolites) == 1

    rxn = ∅ → 1.0nadh
    @test length(rxn.metabolites) == 1
    
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test prod(values(rxn.metabolites)) == -12
    @test ("q8h2_c" in [x.id for x  in keys(rxn.metabolites)])

    rxn = nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0

    @test length(h_p + h_p) == 2
    @test length(h_p + h_p + h_p) == 3

    @test length(+([2.0*q8, nad], nadh)) == 3 
end