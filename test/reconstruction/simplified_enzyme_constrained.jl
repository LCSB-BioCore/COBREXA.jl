@testset "SMOMENT" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # add molar masses to gene products
    for gid in genes(model)
        model.genes[gid].product_molar_mass = get(ecoli_core_gene_product_masses, gid, 0.0)
    end
    model.genes["s0001"] = Gene(id="s0001"; product_molar_mass = 0.0)

    # update isozymes with kinetic information
    for rid in reactions(model)
        if haskey(ecoli_core_reaction_kcats, rid) # if has kcat, then has grr
            newisozymes = Isozyme[]
            for (i, grr) in enumerate(reaction_gene_associations(model, rid))
                push!(
                    newisozymes, 
                    Isozyme(
                        gene_product_stoichiometry = Dict(grr .=> ecoli_core_protein_stoichiometry[rid][i]),
                        kcat_forward = ecoli_core_reaction_kcats[rid][i][1],
                        kcat_backward = ecoli_core_reaction_kcats[rid][i][2]
                    )
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
            lower = [-1000.0, -1.0],
            upper = [nothing, 12.0],
        ) |>
        with_simplified_enzyme_constrained(
            total_enzyme_capacity = 100.0,
        )

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
