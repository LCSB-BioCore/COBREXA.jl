@testset "atom exchanges" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    fluxes = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    # atom tracker
    atom_fluxes = atom_exchange(model, fluxes)
    @test isapprox(atom_fluxes["C"], 0; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes["O"], 0; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes["N"], 0; atol = TEST_TOLERANCE)
end
