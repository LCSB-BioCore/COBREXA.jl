@testset "Constrained allocation FBA" begin

    m = ObjectModel()

    add_reactions!(
        m,
        [
            ReactionForward("r1", Dict("m1" => 1)),
            ReactionForward("r2", Dict("m2" => 1)),
            ReactionForward("r3", Dict("m1" => -1, "m2" => -1, "m3" => 1)),
            ReactionForward("r4", Dict("m3" => -1, "m4" => 1)),
            ReactionBidirectional("r5", Dict("m2" => -1, "m4" => 1)),
            ReactionBidirectional("r6", Dict("m4" => -1)), # different!
        ],
    )

    m.reactions["r3"].gene_associations =
        [Isozyme(["g1"]; kcat_forward = 1.0, kcat_backward = 1.0)]
    m.reactions["r4"].gene_associations = [
        Isozyme(["g1"]; kcat_forward = 2.0, kcat_backward = 2.0),
        Isozyme(["g2"]; kcat_forward = 3.0, kcat_backward = 3.0),
    ]
    m.reactions["r5"].gene_associations = [
        Isozyme(;
            gene_product_stoichiometry = Dict("g3" => 1, "g4" => 2),
            kcat_forward = 70.0,
            kcat_backward = 70.0,
        ),
    ]
    m.reactions["r6"].gene_associations = [ # this should get removed
        Isozyme(;
            gene_product_stoichiometry = Dict("g3" => 1),
            kcat_forward = 10.0,
            kcat_backward = 0.0,
        ),
    ]
    m.objective = Dict("r6" => 1.0)

    add_genes!(
        m,
        [
            Gene(id = "g$i", product_upper_bound = 10.0, product_molar_mass = float(i)) for
            i = 1:4
        ],
    )

    add_metabolites!(m, [Metabolite("m$i") for i = 1:4])

    ribomodel = m |> with_virtualribosome("r6", 0.2)

    @test haskey(ribomodel.genes, "virtualribosome")
    @test first(ribomodel.reactions["r6"].gene_associations).kcat_forward == 0.2
    @test first(m.reactions["r6"].gene_associations).kcat_forward == 10.0


    cam = make_simplified_enzyme_constrained_model(ribomodel; total_enzyme_capacity = 0.5)

    @test coupling(cam)[1, 7] == 5.0

    rxn_fluxes =
        flux_balance_analysis(
            cam,
            Tulip.Optimizer;
            modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
        ) |> values_dict(:reaction)
    @test isapprox(rxn_fluxes["r6"], 0.09695290851008717, atol = TEST_TOLERANCE)

    # test inplace variant
    add_virtualribosome!(m, "r6", 0.2)
    cam = make_simplified_enzyme_constrained_model(m; total_enzyme_capacity = 0.5)

    @test coupling(cam)[1, 7] == 5.0

    rxn_fluxes =
        flux_balance_analysis(
            cam,
            Tulip.Optimizer;
            modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
        ) |> values_dict(:reaction)
    @test isapprox(rxn_fluxes["r6"], 0.09695290851008717, atol = TEST_TOLERANCE)

    # test with_isozyme functions
    iso1 = Isozyme(["g1"]; kcat_forward = 200.0, kcat_backward = 300.0)
    iso2 = Isozyme(["g2"]; kcat_forward = 100.0, kcat_backward = 500.0)
    m2 = m |> with_isozymes(["r3", "r4"], [[iso1], [iso2]])
    @test first(m2.reactions["r3"].gene_associations).kcat_backward == 300.0
    @test first(m2.reactions["r4"].gene_associations).kcat_backward == 500.0
    @test first(m.reactions["r3"].gene_associations).kcat_backward == 1.0

    m2 = m |> with_isozymes("r3", [iso2])
    @test first(m2.reactions["r3"].gene_associations).kcat_backward == 500.0
    @test first(m.reactions["r3"].gene_associations).kcat_backward == 1.0

    add_isozymes!(m, ["r3", "r4"], [[iso1], [iso2]])
    @test first(m.reactions["r3"].gene_associations).kcat_backward == 300.0
    @test first(m.reactions["r4"].gene_associations).kcat_backward == 500.0

end
