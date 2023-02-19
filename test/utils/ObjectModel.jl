@testset "ObjectModel utilities" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # FBA
    fluxes =
        flux_balance_analysis(
            model,
            Tulip.Optimizer;
            modifications = [change_objective("BIOMASS_Ecoli_core_w_GAM")],
        ) |> values_dict

    # bounds setting # TODO why is this here?
    cbm = make_optimization_model(model, Tulip.Optimizer)
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    glucose_index = first(indexin(["EX_glc__D_e"], variables(model)))
    o2_index = first(indexin(["EX_o2_e"], variables(model)))
    atpm_index = first(indexin(["ATPM"], variables(model)))
    set_optmodel_bound!(glucose_index, cbm; upper_bound = -1.0, lower_bound = -1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_optmodel_bound!(o2_index, cbm; upper_bound = 1.0, lower_bound = 1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0
end
