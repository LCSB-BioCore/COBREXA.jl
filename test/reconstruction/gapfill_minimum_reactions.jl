@testset "Gap fill with minimum reactions" begin
    #=
    Implement the small model that should be gapfilled.
    =#
    model = ObjectModel(id = "partial model")

    (m1, m2, m3, m4, m5, m6, m7, m8) = [Metabolite(id = "m$i") for i = 1:8]

    #=
    here reactions are added to the model, but some are commented out. The goal
    of the gap filling is to identify these commented reactions.
    =#
    add_reactions!(
        model,
        [
            ReactionForward("r1", Dict("m1" => 1); default_bound = 1),
            ReactionBidirectional("r2", Dict("m1" => -1, "m2" => 1); default_bound = 10),
            ReactionForward("r3", Dict("m1" => -1, "m3" => 1)),
            ReactionBidirectional("r4", Dict("m2" => -1, "m4" => 1)),
            # ReactionForward("r5", Dict("m3" => -1, "m4" => 1)),
            ReactionForward("r6", Dict("m4" => -1)),
            # ReactionForward("r7", Dict("m2" => -1, "m7" => 1, "m6" => 1)),
            ReactionForward("r8", Dict("m7" => -1, "m8" => 1)),
            ReactionForward("r9", Dict("m8" => -1)),
            # ReactionForward("r10", Dict("m6" => -1)),
            ReactionForward("r11", Dict("m2" => -1, "m3" => -1, "m7" => -1)),
            ReactionBidirectional(
                "r12",
                Dict("m3" => -1, "m5" => 1);
                default_bound = 10,
            ),
        ],
    )

    model.objective = Dict("r11" => 1)

    add_metabolites!(model, [m1, m2, m3, m4, m5, m7, m8])

    r5 = ReactionForward("r5", Dict("m3" => -1, "m4" => 1))
    r7 = ReactionForward("r7", Dict("m2" => -1, "m7" => 1, "m6" => 1))
    r10 = ReactionForward("r10", Dict("m6" => -1))
    rA = ReactionForward("rA", Dict("m1" => -1, "m2" => 1, "m3" => 1))
    rB = ReactionForward("rB", Dict("m2" => -1, "m9" => 1))
    rC = ReactionBidirectional("rC", Dict("m9" => -1, "m10" => 1))
    rD = ReactionBackward("rC", Dict("m10" => -1))

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
