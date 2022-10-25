@testset "Gap fill with minimum reactions" begin
    #=
    Implement the small model that should be gapfilled.
    =#
    model = ObjectModel("partial model")

    (m1, m2, m3, m4, m5, m6, m7, m8) = Metabolite.("m$i" for i = 1:8)

    @add_reactions! model begin
        "r1", nothing → m1, 0, 1
        "r2", m1 ↔ m2, -10, 100
        "r3", m1 → m3, 0, 100
        "r4", m2 ↔ m4, 0, 100
        # "r5", m3 → m4, 0, 100
        "r6", m4 → nothing, 0, 100
        # "r7", m2 → m7 + m6, 0, 100
        "r8", m7 → m8, 0, 100
        "r9", m8 → nothing, 0, 100
        # "r10", m6 → nothing, 0, 100
        "r11", m2 + m3 + m7 → nothing, 0, 100
        "r12", m3 → m5, -10, 10
    end

    model.reactions["r11"].objective_coefficient = 1.0

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
