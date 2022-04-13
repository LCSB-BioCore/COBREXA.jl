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

    get_reaction_isozyme_masses =
        rid ->
            haskey(ecoli_core_protein_stoichiometry, rid) ?
            [
                sum(
                    counts .*
                    get.(Ref(ecoli_core_protein_masses), gids, 0.0),
                ) for (gids, counts) in zip(reaction_gene_association(model, rid), ecoli_core_protein_stoichiometry[rid])
            ] : []

    total_protein_mass = 100.0

    gm =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"];
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        ) |>
        with_gecko(
            reaction_isozymes = get_reaction_isozymes,
            reaction_isozyme_masses = get_reaction_isozyme_masses,
            gene_product_limit = g -> g == "b2779" ? (0.01, 0.06) : (0.0, 1.0),
            mass_fraction_limit = _ -> total_protein_mass,
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
