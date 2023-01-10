@testset "SMOMENT" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])
    for gid in genes(model)
        model.genes[gid].protein_molar_mass = get(ecoli_core_gene_product_masses, gid, 0.0)
    end
    for rid in reactions(model)
        if haskey(ecoli_core_reaction_kcats, rid) # if has kcat, then has grr
            model.reactions[rid].kcat_forward = ecoli_core_reaction_kcats[rid][1][1] # all kcat forwards the same
            model.reactions[rid].kcat_backward = ecoli_core_reaction_kcats[rid][1][2] # all kcat forwards the same
            newisozymes = Isozyme[]
            for (i, grr) in enumerate(reaction_gene_associations(model, rid))
                push!(newisozymes, Isozyme(stoichiometry = Dict(grr .=> ecoli_core_protein_stoichiometry[rid][i])))
            end
            model.reactions[rid].gene_associations = newisozymes
        end
    end

    simplified_enzyme_constrained_model =
        model |>
        with_changed_bounds(
            ["EX_glc__D_e", "GLCpts"],
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        ) |>
        with_simplified_enzyme_constrained(
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
