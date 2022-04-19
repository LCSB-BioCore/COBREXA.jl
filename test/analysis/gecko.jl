@testset "GECKO" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    get_reaction_isozymes =
        rid ->
            haskey(ecoli_core_reaction_kcats, rid) ?
            collect(
                Isozyme(
                    Dict(grr .=> ecoli_core_protein_stoichiometry[rid][i]),
                    ecoli_core_reaction_kcats[rid][i]...,
                ) for (i, grr) in enumerate(reaction_gene_association(model, rid))
            ) : Isozyme[]

    get_gene_product_mass = gid -> get(ecoli_core_protein_masses, gid, 0.0)

    total_protein_mass = 100.0

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
            gene_mass_group_bound = _ -> total_protein_mass,
        )

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = flux_dict(gm, opt_model)
    prot_concens = protein_dict(gm, opt_model)

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.812827846796761,
        atol = TEST_TOLERANCE,
    )

    prot_mass = sum(ecoli_core_protein_masses[gid] * c for (gid, c) in prot_concens)

    @test isapprox(prot_mass, total_protein_mass, atol = TEST_TOLERANCE)
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
        "r1", nothing → m1, -100, 100
        "r2", nothing → m2, -100, 100
        "r3", m1 + m2 → m3, 0, 100
        "r4", m3 ↔ m4, -100, 100 # make reversible instead
        "r5", m2 ↔ m4, -100, 100
        "r6", nothing → m4, 0, 100 
    end

    gs = [Gene("g$i") for i in 1:4]

    m.reactions["r3"].grr = [["g1"]]
    m.reactions["r4"].grr = [["g1"], ["g2"]]
    m.reactions["r5"].grr = [["g3", "g4"]]
    m.reactions["r4"].objective_coefficient = 1.0

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4])

    reaction_isozymes = Dict(
        "r3" => [
            Isozyme(
                Dict("g1" => 1),
                1.0, 
                1.0,
            ),
        ],
        "r4" => [
            Isozyme(
                Dict("g1" => 1),
                2.0, 
                2.0,
            ),
            Isozyme( 
                Dict("g2" => 1),
                3.0, 
                3.0,
            ),
        ],
        "r5" => [
            Isozyme(
                Dict("g3" => 1, "g4" => 2),
                5.0, 
                5.0,
            ),
        ],
    )
    gene_product_bounds = Dict(
        "g1" => (0.0, 10.0),
        "g2" => (0.0, 10.0),
        "g3" => (0.0, 10.0),
        "g4" => (0.0, 10.0),
    )

    gene_product_molar_mass = Dict(
        "g1" => 1.0,
        "g2" => 2.0,
        "g3" => 3.0,
        "g4" => 4.0,
    )

    gene_mass_group_bound = Dict("uncategorized" => 0.5)

    gm =  make_gecko_model(
        m;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_mass_group_bound,
    )

    S = stoichiometry(gm)
    l, u = bounds(gm)

    coupling(gm)
    cl, cu = coupling_bounds(gm)

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = flux_dict(gm, opt_model)
    gene_products = protein_dict(gm, opt_model)
    mass_groups = protein_mass_group_dict(gm, opt_model)

    @test isapprox(rxn_fluxes["r4"], 0.142857, atol = TEST_TOLERANCE)
    @test isapprox(gene_products["g3"], 0.0285714, atol = TEST_TOLERANCE)
    @test isapprox(mass_groups["uncategorized"], 0.5, atol = TEST_TOLERANCE)
end