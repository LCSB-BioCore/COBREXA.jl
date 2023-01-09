@testset "SMOMENT" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    get_gene_product_mass = gid -> get(ecoli_core_gene_product_masses, gid, 0.0)

    get_reaction_isozyme =
        rid ->
            haskey(ecoli_core_reaction_kcats, rid) ?
            argmax(
                simplified_enzyme_constrained_isozyme_speed(get_gene_product_mass),
                Isozyme(
                    stoichiometry = Dict(grr .=> ecoli_core_protein_stoichiometry[rid][i]),
                    kcat_forward = ecoli_core_reaction_kcats[rid][i][1],
                    kcat_backward = ecoli_core_reaction_kcats[rid][i][2],
                ) for (i, grr) in enumerate(reaction_gene_association(model, rid))
            ) : nothing

    simplified_enzyme_constrained_model =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"],
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        ) |>
        with_simplified_enzyme_constrained(
            reaction_isozyme = get_reaction_isozyme,
            gene_product_molar_mass = get_gene_product_mass,
            total_enzyme_capacity = 100.0,
        )
    objective(simplified_enzyme_constrained_model)

    rxn_fluxes = flux_balance_analysis_dict(
        simplified_enzyme_constrained_model,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8907273630431708,
        atol = TEST_TOLERANCE,
    )
end
