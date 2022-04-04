@testset "SMOMENT" begin
    smodel = load_model(StandardModel, model_paths["e_coli_core.json"])
    smodel.reactions["EX_glc__D_e"].lb = -1000.0 # unconstrain because enzyme constraints take over
    flux_measurements = Dict("GLCpts" => (-1.0, 12.0))
    total_protein_mass = 100 # mg/gdW

    remove_slow_isozymes!(
        smodel;
        reaction_protein_stoichiometry = ecoli_core_protein_stoichiometry,
        protein_masses = ecoli_core_protein_masses,
        reaction_kcats = ecoli_core_reaction_kcats,
    )

    model = SMomentModel(
        smodel;
        reaction_kcats = ecoli_core_reaction_kcats,
        reaction_protein_stoichiometry = ecoli_core_protein_stoichiometry,
        protein_masses = ecoli_core_protein_masses,
        total_protein_mass = total_protein_mass, # mg/gdW
        flux_measurements,
    )

    opt_model = flux_balance_analysis(
        model,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
        sense = COBREXA.MOI.MAX_SENSE,
    )

    rxn_fluxes = flux_dict(model, opt_model)

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8907273630431708,
        atol = TEST_TOLERANCE,
    )
end
