@testset "GECKO" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])
    model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstraint because enzyme constraints take over
    total_protein_mass = 100 # mg/gdW

    rxn_fluxes, prot_concens = gecko(
        model,
        Tulip.Optimizer;
        objective_id = "BIOMASS_Ecoli_core_w_GAM§FOR",
        protein_stoichiometry = ecoli_core_protein_stoichiometry,
        protein_masses = ecoli_core_protein_masses,
        reaction_kcats = ecoli_core_reaction_kcats,
        lb_protein_measurements = Dict("b2779" => 0.01),
        ub_protein_measurements = Dict("b2779" => 0.06),
        lb_flux_measurements = Dict("GLCpts" => -1.0),
        ub_flux_measurements = Dict("GLCpts" => 12.0),
        total_protein_mass,
        sense = COBREXA.MOI.MAX_SENSE,
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
        ]
    )

    prot_mass = sum(ecoli_core_protein_masses[gid] * c for (gid, c) in prot_concens)
    @test isapprox(rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"], 0.812827846796761, atol = TEST_TOLERANCE)
    @test isapprox(prot_mass, total_protein_mass, atol = TEST_TOLERANCE)
end



