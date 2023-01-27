@testset "SMOMENT" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    get_gene_product_mass = gid -> get(ecoli_core_gene_product_masses, gid, 0.0)

    get_reaction_isozyme =
        rid ->
            haskey(ecoli_core_reaction_kcats, rid) ?
            argmax(
                smoment_isozyme_speed(get_gene_product_mass),
                Isozyme(
                    Dict(grr .=> fill(1.0, size(grr))),
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
