@testset "GECKO" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # add molar masses to gene products
    for gid in genes(model)
        model.genes[gid].product_molar_mass = get(ecoli_core_gene_product_masses, gid, 0.0)
    end
    model.genes["s0001"] = Gene(id = "s0001"; product_molar_mass = 0.0)

    # update isozymes with kinetic information
    for rid in reactions(model)
        if haskey(ecoli_core_reaction_kcats, rid) # if has kcat, then has grr
            newisozymes = Isozyme[]
            for (i, grr) in enumerate(reaction_gene_associations(model, rid))
                push!(
                    newisozymes,
                    Isozyme(
                        gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))),
                        kcat_forward = ecoli_core_reaction_kcats[rid][i][1],
                        kcat_backward = ecoli_core_reaction_kcats[rid][i][2],
                    ),
                )
            end
            model.reactions[rid].gene_associations = newisozymes
        else
            model.reactions[rid].gene_associations = nothing
        end
    end

    total_enzyme_capacity = 100.0

    # set gene product bounds
    for gid in genes(model)
        lb, ub = gid == "b2779" ? (0.01, 0.06) : (0.0, 1.0)
        model.genes[gid].product_lower_bound = lb
        model.genes[gid].product_upper_bound = ub
    end

    gm =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"];
            lower_bounds = [-1000.0, -1.0],
            upper_bounds = [nothing, 12.0],
        ) |>
        with_enzyme_constraints(; total_enzyme_capacity)

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = values_dict(:reaction, gm, opt_model)
    prot_concens = values_dict(:enzyme, gm, opt_model)

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8128427019072836,
        atol = TEST_TOLERANCE,
    )

    mass_groups = values_dict(:enzyme_group, gm, opt_model)

    @test isapprox(
        sum(values(prot_concens)),
        total_gene_product_mass,
        atol = TEST_TOLERANCE,
    )
    @test isapprox(
        sum(values(prot_concens)),
        mass_groups["uncategorized"],
        atol = TEST_TOLERANCE,
    )

    # test enzyme objective
    growth_lb = rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"] * 0.9
    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [
            change_objective(genes(gm); weights = [], sense = MIN_SENSE),
            change_constraint("BIOMASS_Ecoli_core_w_GAM", lower_bound = growth_lb),
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
        ],
    )
    mass_groups_min = values_dict(:enzyme_group, gm, opt_model)
    @test mass_groups_min["uncategorized"] < mass_groups["uncategorized"]
end

@testset "GECKO small model" begin
    #=
    Implement the small model found in the supplment of the
    original GECKO paper. This model is nice to troubleshoot with,
    because the stoich matrix is small.
    =#
    m = ObjectModel()
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")

    add_reactions!(
        m,
        [
            ReactionForward("r1", Dict("m1" => 1)),
            ReactionForward("r2", Dict("m2" => 1)),
            ReactionForward("r3", Dict("m1" => -1, "m2" => -1, "m3" => 1)),
            ReactionForward("r4", Dict("m3" => -1, "m4" => 1)),
            ReactionBidirectional("r5", Dict("m2" => -1, "m4" => 1)),
            ReactionForward("r6", Dict("m4" => -1)),
        ],
    )

    gs = [
        Gene(id = "g1", product_upper_bound = 10.0, product_molar_mass = 1.0)
        Gene(id = "g2", product_upper_bound = 10.0, product_molar_mass = 2.0)
        Gene(id = "g3", product_upper_bound = 10.0, product_molar_mass = 3.0)
        Gene(id = "g4", product_upper_bound = 10.0, product_molar_mass = 4.0)
    ]

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
    m.objective = Dict("r6" => 1.0)

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4])

    gm = make_enzyme_constrained_model(
        m;
        gene_product_mass_group = Dict("uncategorized" => genes(m), "bound2" => ["g3"]),
        gene_product_mass_group_bound = Dict("uncategorized" => 0.5, "bound2" => 0.04),
    )

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    rxn_fluxes = values_dict(:reaction, gm, opt_model)
    gene_products = values_dict(:enzyme, gm, opt_model)
    mass_groups = values_dict(:enzyme_group, gm, opt_model)

    @test isapprox(rxn_fluxes["r6"], 1.1688888886502442, atol = TEST_TOLERANCE)
    @test isapprox(gene_products["g4"], 0.02666666666304931 * 4.0, atol = TEST_TOLERANCE)
    @test isapprox(mass_groups["uncategorized"], 0.5, atol = TEST_TOLERANCE)
    @test isapprox(mass_groups["bound2"], 0.04, atol = TEST_TOLERANCE)
    @test length(genes(gm)) == 4
    @test length(genes(gm.inner)) == 4
end
