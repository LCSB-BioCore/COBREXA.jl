"""
This following system allows for the thermodynamically infeasible loop
A -> C -> B -> A as the reaction A <-> B is reversible
Reaction                   
---------------------------
∅ -v1-> A                  
        A <-v2-> B         
        A  -v3-> C         
        C  -v4-> B         
                 B -v5-> ∅ 
If we constrain the influx to 1, we can that in the vanilla fba solution by the
fluxes that are > 1, e.g. The exact numerical values don't matter, it's just
important that they are greater than the influx
| Rxn | Flux |
| --- | ---  |
| v1  | 1    |
| v2  | -267 |
| v3  | 268  |
| v4  | 268  |
| v5  | 1    |
"""
function create_loopless_test_model()
    model = CobraModel()
    A = Metabolite("A")
    B = Metabolite("B")
    C = Metabolite("C")
    D = Metabolite("D")

    v1 = ∅ ⟶ A
    v1.id = "v1"
    v1.ub = 1

    v2 = A ⟷ B
    v2.id = "v2"

    v3 = A ⟶ C
    v3.id = "v3"

    v4 = C ⟶ B
    v4.id = "v4"

    v5 = B ⟶ ∅
    v5.id = "v5"

    add!(model, A)
    add!(model, B)
    add!(model, C)
    add!(model, v1)
    add!(model, v2)
    add!(model, v3)
    add!(model, v4)
    add!(model, v5)
    return model
end


@testset "loopless" begin
    model = create_loopless_test_model()
    optimizer = Tulip.Optimizer
    # first test that fba in fact does not give the loopless solution
    objective = findfirst(model.reactions, "v5")
    res_fba = fba(model, optimizer, objective_func=objective)
    # Check if there is any flux that exceeds the input
    @test !all(abs.(values(res_fba)) .<= 1)

    # Now check that the flux we get the right solution
    res_loopless = loopless(model, optimizer, objective)
    # Check that there is no flux that exceeds the input
    @test all(abs.(values(res_loopless)) .<= 1)
end
