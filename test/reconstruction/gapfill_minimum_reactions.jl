@testset "Gap fill with minimum reactions" begin
    #=
    Implement the small model that should be gapfilled.
    =#
    model = ObjectModel(id = "partial model")

    (m1, m2, m3, m4, m5, m6, m7, m8) = [Metabolite(id = "m$i") for i = 1:8]

    add_reactions!(
        model,
        [
            # Reaction("r1", Dict("m1" => 1), :forward; upper_bound = 1),
            Reaction("r2", Dict("m1" => -1, "m2" => 1), :bidirectional; lower_bound = -10),
            Reaction("r3", Dict("m1" => -1, "m3" => 1), :forward),
            Reaction("r4", Dict("m2" => -1, "m4" => 1), :bidirectional),
            # Reaction("r5", Dict("m3" => -1, "m4" => 1), :forward),
            Reaction("r6", Dict("m4" => -1), :forward),
            # Reaction("r7", Dict("m2" => -1, "m7" => 1, "m6" => 1), :forward),
            Reaction("r8", Dict("m7" => -1, "m8" => 1), :forward),
            Reaction("r9", Dict("m8" => -1), :forward),
            # Reaction("r10", Dict("m6" => -1), :forward),
            Reaction("r11", Dict("m2" => -1, "m3" => -1, "m7" => -1), :forward),
            Reaction("r12", Dict("m3" => -1, "m5" => 1), :forward; lower_bound = -10, upper_bound = 10),
            
        ]
    )

    model.objective = Dict("r11" => 1)

    add_metabolites!(model, [m1, m2, m3, m4, m5, m7, m8])

    r5 = Reaction("r5", Dict("m3" => -1, "m4" => 1), :forward)
    r7 = Reaction("r7", Dict("m2" => -1, "m7" => 1, "m6" => 1), :forward)
    r10 = Reaction("r10", Dict("m6" => -1), :forward)
    rA = Reaction("rA", Dict("m1" => -1, "m2" => 1, "m3" => 1), :forward)
    rB = Reaction("rB", Dict("m2" => -1, "m9" => 1), :forward)
    rC = Reaction("rC", Dict("m9" => -1, "m10" => 1), :bidirectional)
    rD = Reaction("rC", Dict("m10" => -1), :reverse)

    universal_reactions = [r5, r7, r10, rA, rB, rC, rD]

    rxns =
        gapfill_minimum_reactions(
            model,
            universal_reactions,
            GLPK.Optimizer;
            objective_bounds = (0.1, 1000.0),
        ) |> gapfilled_rids(universal_reactions)

    @test issetequal(["r7", "r10"], rxns)
end
