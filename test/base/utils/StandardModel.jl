@testset "StandardModel utilities" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    # FBA
    glc = model.reactions["EX_glc__D_e"]
    optimizer = Tulip.Optimizer # quiet by default
    sol = flux_balance_analysis_dict(
        model,
        optimizer;
        modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
    )

    # atom tracker
    atom_fluxes = atom_exchange(sol, model)
    @test isapprox(atom_fluxes["C"], 37.19016648975907; atol = TEST_TOLERANCE)
    @test atom_exchange("FBA", model)["C"] == 0.0
    @test isapprox(
        atom_exchange("BIOMASS_Ecoli_core_w_GAM", model)["C"],
        -42.5555;
        atol = TEST_TOLERANCE,
    )

    # metabolite trackers
    consuming, producing = metabolite_fluxes(sol, model)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738; atol = TEST_TOLERANCE)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818; atol = TEST_TOLERANCE)

    # set bounds
    cbm = make_optimization_model(model, optimizer)
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    glucose_index = first(indexin(["EX_glc__D_e"], reactions(model)))
    o2_index = first(indexin(["EX_o2_e"], reactions(model)))
    atpm_index = first(indexin(["ATPM"], reactions(model)))
    set_optmodel_bound!(glucose_index, cbm; ub = -1.0, lb = -1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_optmodel_bound!(o2_index, cbm; ub = 1.0, lb = 1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0
end
