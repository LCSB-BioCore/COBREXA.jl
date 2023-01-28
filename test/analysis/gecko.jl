@testset "GECKO" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    get_reaction_isozymes =
        rid ->
            haskey(ecoli_core_reaction_kcats, rid) ?
            collect(
                Isozyme(
                    Dict(grr .=> fill(1.0, size(grr))),
                    ecoli_core_reaction_kcats[rid][i]...,
                ) for (i, grr) in enumerate(reaction_gene_association(model, rid))
            ) : Isozyme[]

    get_gene_product_mass = gid -> get(ecoli_core_gene_product_masses, gid, 0.0)

    total_gene_product_mass = 100.0

    bounded_model =
        model |> with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"];
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        )

    gm =
        bounded_model |> with_gecko(
            reaction_isozymes = get_reaction_isozymes,
            gene_product_bounds = g -> g == "b2779" ? (0.01, 0.06) : (0.0, 1.0),
            gene_product_molar_mass = get_gene_product_mass,
            gene_product_mass_group_bound = _ -> total_gene_product_mass,
        )

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = flux_dict(gm, opt_model)
    prot_concens = gene_product_dict(gm, opt_model)

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8129179015245396,
        atol = TEST_TOLERANCE,
    )

    prot_mass = sum(ecoli_core_gene_product_masses[gid] * c for (gid, c) in prot_concens)
    mass_groups = gene_product_mass_group_dict(gm, opt_model)

    @test isapprox(prot_mass, total_gene_product_mass, atol = TEST_TOLERANCE)
    @test isapprox(prot_mass, mass_groups["uncategorized"], atol = TEST_TOLERANCE)

    # test enzyme objective
    growth_lb = rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"] * 0.9
    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [
            change_objective(genes(gm); weights = [], sense = COBREXA.MIN_SENSE),
            change_constraint("BIOMASS_Ecoli_core_w_GAM", lb = growth_lb),
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
        ],
    )
    mass_groups_min = gene_product_mass_group_dict(gm, opt_model)
    @test mass_groups_min["uncategorized"] < mass_groups["uncategorized"]
end

@testset "GECKO small model" begin
    #=
    Implement the small model found in the supplment of the
    original GECKO paper. This model is nice to troubleshoot with,
    because the stoich matrix is small.
    =#
    m = StandardModel("gecko")
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")

    @add_reactions! m begin
        "r1", nothing → m1, 0, 100
        "r2", nothing → m2, 0, 100
        "r3", m1 + m2 → m3, 0, 100
        "r4", m3 → m4, 0, 100
        "r5", m2 ↔ m4, -100, 100
        "r6", m4 → nothing, 0, 100
    end

    gs = [Gene("g$i") for i = 1:5]

    m.reactions["r2"].grr = [["g5"]]
    m.reactions["r3"].grr = [["g1"]]
    m.reactions["r4"].grr = [["g1"], ["g2"]]
    m.reactions["r5"].grr = [["g3", "g4"]]
    m.reactions["r6"].objective_coefficient = 1.0

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4])

    reaction_isozymes = Dict(
        "r3" => [Isozyme(Dict("g1" => 1), 1.0, 1.0)],
        "r4" =>
            [Isozyme(Dict("g1" => 1), 2.0, 2.0), Isozyme(Dict("g2" => 1), 3.0, 3.0)],
        "r5" => [Isozyme(Dict("g3" => 1, "g4" => 2), 70.0, 70.0)],
    )
    gene_product_bounds = Dict(
        "g1" => (0.0, 10.0),
        "g2" => (0.0, 10.0),
        "g3" => (0.0, 10.0),
        "g4" => (0.0, 10.0),
    )

    gene_product_molar_mass = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0)

    gene_product_mass_group_bound = Dict("uncategorized" => 0.5)

    gm = make_gecko_model(
        m;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = flux_dict(gm, opt_model)
    gene_products = gene_product_dict(gm, opt_model)
    mass_groups = gene_product_mass_group_dict(gm, opt_model)

    @test isapprox(rxn_fluxes["r6"], 3.181818181753438, atol = TEST_TOLERANCE)
    @test isapprox(gene_products["g4"], 0.09090909090607537, atol = TEST_TOLERANCE)
    @test isapprox(mass_groups["uncategorized"], 0.5, atol = TEST_TOLERANCE)
    @test length(genes(gm)) == 4
    @test length(genes(gm.inner)) == 5
end
