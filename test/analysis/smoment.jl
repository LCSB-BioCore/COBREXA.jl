@testset "SMOMENT" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    get_gene_product_mass = gid -> get(ecoli_core_gene_product_masses, gid, 0.0)

    # fix ordering
    model.reactions["ATPS4r"].grr = [
        ["b3736", "b3737", "b3738", "b3731", "b3732", "b3733", "b3734", "b3735"],
        ["b3736", "b3737", "b3738", "b3731", "b3732", "b3733", "b3734", "b3735", "b3739"]
    ]

    model.reactions["GLCpts"].grr = [
        ["b2417", "b1621", "b2415", "b2416"],
        ["b1817", "b1818", "b1819", "b2415", "b2416"],
        ["b2417", "b1101", "b2415", "b2416"],
    ]

    get_reaction_isozyme =
        rid ->
            haskey(ecoli_core_reaction_kcats, rid) ?
            argmax(
                smoment_isozyme_speed(get_gene_product_mass),
                Isozyme(
                    Dict(grr .=> ecoli_core_protein_stoichiometry[rid][i]),
                    ecoli_core_reaction_kcats[rid][i]...,
                ) for (i, grr) in enumerate(reaction_gene_association(model, rid))
            ) : nothing

    smoment_model =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"],
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        ) |>
        with_smoment(
            reaction_isozyme = get_reaction_isozyme,
            gene_product_molar_mass = get_gene_product_mass,
            total_enzyme_capacity = 100.0,
        )

    rxn_fluxes = flux_balance_analysis_dict(
        smoment_model,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8907347602586123,
        atol = TEST_TOLERANCE,
    )
end
