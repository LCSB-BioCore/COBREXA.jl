
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
    @test isapprox(atom_fluxes["C"], -37.1902, atol = 1e-3)

    # metabolite trackers
    consuming, producing = metabolite_fluxes(sol, model)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738, atol = 1e-3)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818, atol = 1e-3)

    # set bounds
    cbm = make_optimization_model(model, optimizer)
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    glucose_index = first(indexin(["EX_glc__D_e"], reactions(model)))
    o2_index = first(indexin(["EX_o2_e"], reactions(model)))
    atpm_index = first(indexin(["ATPM"], reactions(model)))
    set_bound(glucose_index, cbm; ub = -1.0, lb = -1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_bound(o2_index, cbm; ub = 1.0, lb = 1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0

    # find exchange reactions
    ex_rxns = find_exchange_reactions(model)
    @test length(ex_rxns) == 21
    @test "BIOMASS_Ecoli_core_w_GAM" in ex_rxns
    
    ex_rxns = find_exchange_reactions(model; exclude_biomass=true)
    @test length(ex_rxns) == 20
    @test !("BIOMASS_Ecoli_core_w_GAM" in ex_rxns)
    
    ex_rxn_mets = find_exchange_metabolites(model)
    @test length(ex_rxn_mets) == 21
    @test ex_rxn_mets["EX_for_e"]["for_e"] == -1.0
end
