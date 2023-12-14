
@testset "GECKO small model" begin
    #=
    Implement the small model found in the supplment of the
    original GECKO paper. This model is nice to troubleshoot with,
    because the stoich matrix is small.
    =#
    mets = Dict(
        "m1" => A.CanonicalModel.Metabolite(),
        "m2" => A.CanonicalModel.Metabolite(),
        "m3" => A.CanonicalModel.Metabolite(),
        "m4" => A.CanonicalModel.Metabolite(),
    )

    rxns = Dict(
        "r1" => A.CanonicalModel.Reaction(
            lower_bound = 0.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m1" => 1.0),
        ),
        "r2" => A.CanonicalModel.Reaction(
            lower_bound = 0.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m2" => 1.0),
        ),
        "r3" => A.CanonicalModel.Reaction(
            lower_bound = 0.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => 1.0),
        ),
        "r4" => A.CanonicalModel.Reaction(
            lower_bound = 0.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m3" => -1.0, "m4" => 1.0),
        ),
        "r5" => A.CanonicalModel.Reaction(
            lower_bound = -100.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m2" => -1.0, "m4" => 1.0),
        ),
        "r6" => A.CanonicalModel.Reaction(
            lower_bound = 0.0,
            upper_bound = 100.0,
            stoichiometry = Dict("m4" => -1.0),
            objective_coefficient = 1.0,
        ),
    )

    gs = Dict("g$i" => A.CanonicalModel.Gene() for i = 1:5)

    model = A.CanonicalModel.Model(rxns, mets, gs)

    reaction_isozymes = Dict(
        "r3" => Dict("iso1" => Isozyme(Dict("g1" => 1), 1.0, 1.0)),
        "r4" => Dict(
            "iso1" => Isozyme(Dict("g1" => 1), 2.0, 2.0),
            "iso2" => Isozyme(Dict("g2" => 1), 3.0, 3.0),
        ),
        "r5" => Dict("iso1" => Isozyme(Dict("g3" => 1, "g4" => 2), 70.0, 70.0)),
    )

    gene_product_molar_mass =
        Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)

    m = build_enzyme_constrained_model(
        model,
        reaction_isozymes,
        gene_product_molar_mass,
        [("total_proteome_bound", A.genes(model), 0.5)],
    )

    řešení = optimized_constraints(
        m;
        objective = m.objective.value,
        optimizer = Tulip.Optimizer,
        modifications = [set_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(řešení.objective, 3.181818181753438, atol = TEST_TOLERANCE)
    @test isapprox(řešení.enzymes.g4, 0.09090909090607537, atol = TEST_TOLERANCE)
    @test isapprox(řešení.total_proteome_bound, 0.5, atol = TEST_TOLERANCE)
end
