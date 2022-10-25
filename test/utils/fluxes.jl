@testset "Flux utilities" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    fluxes = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    consuming, producing = metabolite_fluxes(model, fluxes)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738; atol = TEST_TOLERANCE)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818; atol = TEST_TOLERANCE)

    # remove the biomass production from the fluxes, so that there's some atom
    # disbalance that can be measured
    delete!(fluxes, "BIOMASS_Ecoli_core_w_GAM")

    # atom tracker
    atom_fluxes_out = atom_fluxes(model, fluxes)
    @test isapprox(atom_fluxes_out["C"], 37.190166489763214; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes_out["O"], 41.663071522672226; atol = TEST_TOLERANCE)
    @test isapprox(atom_fluxes_out["N"], 4.765319167566247; atol = TEST_TOLERANCE)
end
