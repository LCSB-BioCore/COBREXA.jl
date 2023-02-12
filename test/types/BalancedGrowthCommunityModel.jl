@testset "BalancedGrowthCommunityModel: simple model" begin
    m1 = ObjectModel(id = "Model1")
    add_metabolites!(
        m1,
        [
            Metabolite("A"),
            Metabolite("B"),
            Metabolite("Ae"),
            Metabolite("Be"),
            Metabolite("X1"),
        ],
    )
    add_genes!(m1, [Gene("g1"), Gene("g2"), Gene("g3"), Gene("g4")])
    add_reactions!(
        m1,
        [
            ReactionBidirectional("EX_A", Dict("Ae" => -1)),
            ReactionBidirectional("r1", Dict("Ae" => -1, "A" => 1)),
            ReactionBidirectional("r2", Dict("A" => -1, "B" => 1, "X1" => 1)),
            ReactionForward("r3", Dict("B" => -1, "Be" => 1)),
            ReactionForward("EX_B", Dict("Be" => -1)),
        ],
    )

    m2 = ObjectModel(id = "Model2")
    add_metabolites!(
        m2,
        [
            Metabolite("Ae"),
            Metabolite("A"),
            Metabolite("C"),
            Metabolite("Ce"),
            Metabolite("X2"),
        ],
    )
    add_genes!(m2, [Gene("g1"), Gene("g2"), Gene("g3"), Gene("g4")])
    add_reactions!(
        m2,
        [
            ReactionForward("r3", Dict("C" => -1, "Ce" => 1)),
            ReactionForward("EX_C", Dict("Ce" => -1)),
            ReactionBidirectional("EX_A", Dict("Ae" => -1)),
            ReactionBidirectional("r1", Dict("Ae" => -1, "A" => 1)),
            ReactionBidirectional("r2", Dict("A" => -1, "C" => 1, "X2" => 1)),
        ],
    )

    cm1 = CommunityMember(
        id = "m1",
        abundance = 0.2,
        model = m1,
        exchange_reaction_ids = ["EX_A", "EX_B"],
        biomass_metabolite_id = "X1",
    )
    @test contains(sprint(show, MIME("text/plain"), cm1), "community member")

    cm2 = CommunityMember(
        id = "m2",
        abundance = 0.8,
        model = m2,
        exchange_reaction_ids = ["EX_A", "EX_C"],
        biomass_metabolite_id = "X2",
    )

    cm = BalancedGrowthCommunityModel(
        members = [cm1, cm2],
        env_met_flux_bounds = Dict("Ae" => (-10, 10)),
    )
    @test contains(sprint(show, MIME("text/plain"), cm), "balanced growth")

    @test issetequal(
        variables(cm),
        [
            "m1#EX_A"
            "m1#r1"
            "m1#r2"
            "m1#r3"
            "m1#EX_B"
            "m2#r3"
            "m2#EX_C"
            "m2#EX_A"
            "m2#r1"
            "m2#r2"
            "EX_Ae"
            "EX_Be"
            "EX_Ce"
            "equal_growth_rates_biomass_function"
        ],
    )

    @test issetequal(
        metabolites(cm),
        [
            "m1#A"
            "m1#B"
            "m1#Ae"
            "m1#Be"
            "m1#X1"
            "m2#Ae"
            "m2#A"
            "m2#C"
            "m2#Ce"
            "m2#X2"
            "ENV_Ae"
            "ENV_Be"
            "ENV_Ce"
        ],
    )

    @test issetequal(
        genes(cm),
        [
            "m1#g1"
            "m1#g2"
            "m1#g3"
            "m1#g4"
            "m2#g1"
            "m2#g2"
            "m2#g3"
            "m2#g4"
        ],
    )

    @test n_variables(cm) == 14
    @test n_metabolites(cm) == 13
    @test n_genes(cm) == 8

    @test all(
        stoichiometry(cm) .== [
            0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            -1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.2
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 -1.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 -0.8
            1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 -1.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0
        ],
    )

    lbs, ubs = bounds(cm)
    @test all(
        lbs .== [
            -200.0
            -200.0
            -200.0
            0.0
            0.0
            0.0
            0.0
            -800.0
            -800.0
            -800.0
            -10.0
            -1000.0
            -1000.0
            0.0
        ],
    )
    @test all(
        ubs .== [
            200.0
            200.0
            200.0
            200.0
            200.0
            800.0
            800.0
            800.0
            800.0
            800.0
            10.0
            1000.0
            1000.0
            1000.0
        ],
    )

    @test all(objective(cm) .== [
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        1.0
    ])

    @test n_coupling_constraints(cm) == 0
    @test isempty(coupling(cm))
    @test all(isempty.(coupling_bounds(cm)))
end

@testset "BalancedGrowthCommunityModel: e coli core" begin
    ecoli = load_model(ObjectModel, model_paths["e_coli_core.json"])

    add_biomass_metabolite!(ecoli, "BIOMASS_Ecoli_core_w_GAM")

    ecoli.reactions["EX_glc__D_e"].lower_bound = -1000

    ex_rxns = COBREXA.Utils.find_exchange_reaction_ids(ecoli)

    a1 = 0.2 # abundance species 1
    a2 = 0.8 # abundance species 2

    cm1 = CommunityMember(
        id = "ecoli1",
        abundance = a1,
        model = ecoli,
        exchange_reaction_ids = ex_rxns,
        biomass_metabolite_id = "biomass",
    )
    cm2 = CommunityMember(
        id = "ecoli2",
        abundance = a2,
        model = ecoli,
        exchange_reaction_ids = ex_rxns,
        biomass_metabolite_id = "biomass",
    )

    cm = BalancedGrowthCommunityModel(
        members = [cm1, cm2],
        env_met_flux_bounds = Dict("glc__D_e" => (-10, 10)),
    )

    d = flux_balance_analysis_dict(cm, Tulip.Optimizer)

    @test isapprox(d[cm.objective_id], 0.8739215069521299, atol = TEST_TOLERANCE)

    @test isapprox(
        d["ecoli1#BIOMASS_Ecoli_core_w_GAM"],
        a1 * d[cm.objective_id],
        atol = TEST_TOLERANCE,
    )

    @test isapprox(
        d["ecoli2#BIOMASS_Ecoli_core_w_GAM"],
        a2 * d[cm.objective_id],
        atol = TEST_TOLERANCE,
    )

    @test isapprox(
        d["EX_glc__D_e"],
        d["ecoli1#EX_glc__D_e"] + d["ecoli2#EX_glc__D_e"],
        atol = TEST_TOLERANCE,
    )

    # test if model can be converted to another type
    om = convert(ObjectModel, cm)
    @test n_variables(om) == n_variables(cm)
end

@testset "BalancedGrowthCommunityModel: enzyme constrained e coli" begin
    ecoli = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # test added biomass metabolite
    modded_ecoli = ecoli |> with_added_biomass_metabolite("BIOMASS_Ecoli_core_w_GAM")
    @test "biomass" in metabolites(modded_ecoli)
    @test !("biomass" in metabolites(ecoli))
    @test haskey(modded_ecoli.reactions["BIOMASS_Ecoli_core_w_GAM"].metabolites, "biomass")
    @test !haskey(ecoli.reactions["BIOMASS_Ecoli_core_w_GAM"].metabolites, "biomass")

    # add molar masses to gene products
    for gid in genes(ecoli)
        ecoli.genes[gid].product_molar_mass = get(ecoli_core_gene_product_masses, gid, 0.0)
        ecoli.genes[gid].product_upper_bound = 10.0
    end
    ecoli.genes["s0001"] = Gene(id = "s0001"; product_molar_mass = 0.0)
    ecoli.genes["s0001"].product_upper_bound = 10.0

    # update isozymes with kinetic information
    for rid in reactions(ecoli)
        if haskey(ecoli_core_reaction_kcats, rid) # if has kcat, then has grr
            newisozymes = Isozyme[]
            for (i, grr) in enumerate(reaction_gene_associations(ecoli, rid))
                push!(
                    newisozymes,
                    Isozyme(
                        gene_product_stoichiometry = Dict(
                            grr .=> ecoli_core_protein_stoichiometry[rid][i],
                        ),
                        kcat_forward = ecoli_core_reaction_kcats[rid][i][1],
                        kcat_backward = ecoli_core_reaction_kcats[rid][i][2],
                    ),
                )
            end
            ecoli.reactions[rid].gene_associations = newisozymes
        else
            ecoli.reactions[rid].gene_associations = nothing
        end
    end

    gm =
        ecoli |>
        with_added_biomass_metabolite("BIOMASS_Ecoli_core_w_GAM") |>
        with_changed_bounds(["EX_glc__D_e"]; lower_bound = [-1000.0], upper_bound = [0]) |>
        with_enzyme_constraints(
            gene_product_mass_group_bound = Dict("uncategorized" => 100.0),
        )

    ex_rxns = find_exchange_reaction_ids(ecoli)

    a1 = 0.2 # abundance species 1
    a2 = 0.8 # abundance species 2

    cm1 = CommunityMember(
        id = "ecoli1",
        abundance = a1,
        model = gm,
        exchange_reaction_ids = ex_rxns,
        biomass_metabolite_id = "biomass",
    )
    cm2 = CommunityMember(
        id = "ecoli2",
        abundance = a2,
        model = gm,
        exchange_reaction_ids = ex_rxns,
        biomass_metabolite_id = "biomass",
    )

    cm = BalancedGrowthCommunityModel(members = [cm1, cm2])

    opt_model = flux_balance_analysis(
        cm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    f_d = flux_dict(cm, opt_model)
    @test isapprox(f_d[cm.objective_id], 0.9210848582802592, atol = TEST_TOLERANCE)
end
