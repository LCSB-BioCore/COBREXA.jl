@testset "Reaction utilities" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # FBA
    fluxes =
        flux_balance_analysis(
            model,
            Tulip.Optimizer;
            modifications = [modify_objective("BIOMASS_Ecoli_core_w_GAM")],
        ) |> values_dict

    # test if reaction is balanced
    @test reaction_mass_balanced(model, "PFL")
    @test !reaction_mass_balanced(model, "BIOMASS_Ecoli_core_w_GAM")
    @test reaction_mass_balanced(model, model.reactions["PFL"])
    @test !reaction_mass_balanced(model, model.reactions["BIOMASS_Ecoli_core_w_GAM"])
    @test reaction_mass_balanced(model, Dict("h_c" => -1.0, "h_e" => 1.0))
    @test !reaction_mass_balanced(model, Dict("h_c" => -1.0, "h2o_c" => 1.0))

    # test if a reaction is a boundary reaction
    @test !is_boundary(model.reactions["PFL"])
    @test is_boundary(model.reactions["EX_glc__D_e"])
    @test !is_boundary(model, "PFL")
    @test is_boundary(model, "EX_glc__D_e")
    @test !is_boundary(model, model.reactions["PFL"])
    @test is_boundary(model, model.reactions["EX_glc__D_e"])

    # single-reaction atom balance
    @test reaction_atom_balance(model, "FBA")["C"] == 0.0
    @test isapprox(
        reaction_atom_balance(model, "BIOMASS_Ecoli_core_w_GAM")["C"],
        -42.5555;
        atol = TEST_TOLERANCE,
    )
    @test reaction_atom_balance(model, model.reactions["FBA"])["C"] == 0.0
    @test isapprox(
        reaction_atom_balance(model, model.reactions["BIOMASS_Ecoli_core_w_GAM"])["C"],
        -42.5555;
        atol = TEST_TOLERANCE,
    )
    @test reaction_atom_balance(model, Dict("h_c" => -1.0, "h2o_c" => 1.0))["H"] == 1.0

    # test if reaction equation can be built back into a sensible reaction string
    req = Dict("coa_c" => -1, "for_c" => 1, "accoa_c" => 1, "pyr_c" => -1)
    rstr_out = stoichiometry_string(req)
    @test occursin("coa_c", split(rstr_out, " = ")[1])
    @test occursin("for", split(rstr_out, " = ")[2])
end
