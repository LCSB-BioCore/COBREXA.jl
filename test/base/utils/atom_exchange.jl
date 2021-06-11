@testset "atom exchanges" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    fluxes = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    # remove the biomass production from the fluxes, so that there's some atom
    # disbalance that can be measured
    delete!(fluxes, "BIOMASS_Ecoli_core_w_GAM")

    # atom tracker
    atom_fluxes = atom_exchange(model, fluxes)
    @test isapprox(atom_fluxes["C"], 37.190166489763214; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes["O"], 41.663071522672226; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes["N"], 4.765319167566247; atol = TEST_TOLERANCE)
end
