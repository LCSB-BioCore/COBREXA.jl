@testset "Max-min driving force analysis" begin

    # This is a relatively standard model of glycolysis.
    # If editing, keep this the same as the corresponding notebook.
    mets = [
        "13dpg",
        "2pg",
        "3pg",
        "adp",
        "atp",
        "dhap",
        "f6p",
        "fdp",
        "g3p",
        "g6p",
        "glc__D",
        "h",
        "h2o",
        "lac__D",
        "nad",
        "nadh",
        "pep",
        "pi",
        "pyr",
    ]

    rxns = Dict(
        "ENO" => Dict("2pg" => -1.0, "h2o" => 1.0, "pep" => 1),
        "FBA" => Dict("fdp" => -1.0, "dhap" => 1.0, "g3p" => 1.0),
        "GAPD" => Dict(
            "g3p" => -1.0,
            "nad" => -1.0,
            "pi" => -1.0,
            "h" => 1.0,
            "nadh" => 1.0,
            "13dpg" => 1.0,
        ),
        "HEX" => Dict(
            "atp" => -1.0,
            "glc__D" => -1.0,
            "g6p" => 1.0,
            "adp" => 1.0,
            "h" => 1.0,
        ),
        "LDH" => Dict(
            "pyr" => -1.0,
            "nadh" => -1.0,
            "h" => -1.0,
            "nad" => 1.0,
            "lac__D" => 1.0,
        ),
        "PFK" =>
            Dict("f6p" => -1.0, "atp" => -1.0, "adp" => 1.0, "h" => 1.0, "fdp" => 1.0),
        "PGI" => Dict("g6p" => -1.0, "f6p" => 1.0),
        "PGK" => Dict("13dpg" => -1.0, "adp" => -1.0, "atp" => 1.0, "3pg" => 1.0),
        "PGM" => Dict("3pg" => -1.0, "2pg" => 1),
        "PYK" =>
            Dict("pep" => -1.0, "adp" => -1.0, "h" => -1.0, "atp" => 1.0, "pyr" => 1.0),
        "TPI" => Dict("dhap" => -1.0, "g3p" => 1.0),
    )

    model = StandardModel("Glycolysis")

    add_metabolites!(model, Metabolite.(mets))
    add_reactions!(
        model,
        collect(Reaction(rid; metabolites = mets) for (rid, mets) in rxns),
    )

    gibbs_energies = Dict(
        "ENO" => -3.81089,
        "FBA" => 22.3932,
        "GAPD" => 4.60271,
        "HEX" => -17.90,
        "LDH" => -23.6803,
        "PFK" => -14.5988,
        "PGI" => 2.6617,
        "PGK" => -19.32,
        "PGM" => -4.52041,
        "PYK" => -27.5833,
        "TPI" => 5.57535,
    )

    res = max_min_driving_force(
        model,
        gibbs_energies,
        Tulip.Optimizer;
        ignore_metabolites = ["h", "h2o"],
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 500)],
        concentration_ratios = Dict(("atp", "adp") => 10.0, ("nadh", "nad") => 0.1),
        constant_concentrations = Dict("pi" => 10e-3),
        concentration_lb = 1e-6,
        concentration_ub = 10e-3,
    )

    expected_energies = Dict(
        "PFK" => -7.870772134961782,
        "HEX" => -12.989138845164854,
        "PGI" => -7.243212703456118,
        "PGK" => -4.890690058164602,
        "LDH" => -15.276010609005,
        "TPI" => -2.1225264622447533,
        "ENO" => -7.0742938630297125,
        "PYK" => -10.796220177254469,
        "FBA" => -2.1225265140886544,
        "GAPD" => -2.122526375820849,
        "PGM" => -6.024520569096858,
    )
    expected_concentrations = Dict(
        "13dpg" => 1.0000000348183199e-6,
        "dhap" => 0.00336372174320087,
        "2pg" => 1.838217836470006e-5,
        "pep" => 4.928030130227643e-6,
        "3pg" => 3.372141166159583e-5,
        "nad" => 0.0012492046634450203,
        "nadh" => 0.000124920466441468,
        "pi" => 0.010000000000000004,
        "fdp" => 0.009999999541658144,
        "pyr" => 0.00043017235652799234,
        "atp" => 0.0012890806409896136,
        "g6p" => 0.0036021859170920294,
        "adp" => 0.00012890806403201158,
        "g3p" => 0.00015073373256424484,
        "f6p" => 6.62674858888359e-5,
        "glc__D" => 4.9684448349558665e-5,
        "lac__D" => 0.0012764690773165621,
    )
    @test isapprox(res.mmdf, 2.122526369934736, atol = TEST_TOLERANCE)
    @test issetequal(keys(res.energies), keys(expected_energies))
    @test issetequal(keys(res.concentrations), keys(expected_concentrations))
    @test all(
        isapprox(res.energies[i], expected_energies[i], atol = TEST_TOLERANCE) for
        i in keys(expected_energies)
    )
    @test all(
        isapprox(res.concentrations[i], expected_concentrations[i], atol = TEST_TOLERANCE)
        for i in keys(expected_concentrations)
    )
end
