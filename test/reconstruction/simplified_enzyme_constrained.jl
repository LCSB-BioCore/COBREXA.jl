@testset "SMOMENT" begin
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

    simplified_enzyme_constrained_model =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"],
            lower_bounds = [-1000.0, -1.0],
            upper_bounds = [nothing, 12.0],
        ) |>
        with_simplified_enzyme_constraints(total_reaction_mass_bound = 100.0)

    rxn_fluxes =
        flux_balance_analysis(
            simplified_enzyme_constrained_model,
            Tulip.Optimizer;
            modifications = [modify_optimizer_attribute("IPM_IterationsLimit", 1000)],
        ) |> values_dict

    @test isapprox(rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"], 0.890731, atol = TEST_TOLERANCE)
end

@testset "Small SMOMENT" begin

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
    m.reactions["r4"].gene_associations =
        [Isozyme(["g1"]; kcat_forward = 2.0, kcat_backward = 2.0)]
    m.reactions["r5"].gene_associations = [
        Isozyme(;
            gene_product_stoichiometry = Dict("g3" => 1, "g4" => 2),
            kcat_forward = 2.0,
            kcat_backward = 2.0,
        ),
    ]
    m.objective = Dict("r6" => 1.0)

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4])

    sm = make_simplified_enzyme_constrained_model(
        m;
        reaction_mass_groups = Dict("b1" => ["r3", "r4"], "b2" => ["r5", "r4"]),
        reaction_mass_group_bounds = Dict("b2" => 0.5, "b1" => 0.2),
    )

    stoichiometry(sm)

    cpling = sum(coupling(sm), dims = 2)
    @test 1.5 in cpling && 11.5 in cpling

    lbs, ubs = coupling_bounds(sm)
    @test all(lbs .== 0.0)
    @test 0.5 in ubs && 0.2 in ubs

    res = flux_balance_analysis(
        sm,
        Tulip.Optimizer;
        modifications = [modify_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(solved_objective_value(res), 0.21212120975836252; atol = TEST_TOLERANCE)
end
