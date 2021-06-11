@testset "atom exchanges" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    fluxes = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    # atom tracker
    atom_fluxes = atom_exchange(fluxes, model)
    @test isapprox(atom_fluxes["C"], 37.19016648975907; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes["N"], 37.19016648975907; atol = TEST_TOLERANCE)
end
