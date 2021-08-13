@testset "Max min driving force analysis" begin

    model = StandardModel("Glycolysis")

    mets = [
        Metabolite("glc__D"),
        Metabolite("g6p"),
        Metabolite("f6p"),
        Metabolite("fdp"),
        Metabolite("dhap"),
        Metabolite("g3p"),
        Metabolite("13dpg"),
        Metabolite("3pg"),
        Metabolite("2pg"),
        Metabolite("pep"),
        Metabolite("pyr"),
        Metabolite("lac__D"),
        Metabolite("nadh"),
        Metabolite("nad"),
        Metabolite("h"),
        Metabolite("atp"),
        Metabolite("adp"),
        Metabolite("h2o"),
        Metabolite("pi"),
    ]

    rxns = [
        Reaction(
            "HEX";
            metabolites = Dict(
                "atp" => -1.0,
                "glc__D" => -1.0,
                "g6p" => 1.0,
                "adp" => 1.0,
                "h" => 1.0,
            ),
        ),
        Reaction("PGI"; metabolites = Dict("g6p" => -1.0, "f6p" => 1.0)),
        Reaction(
            "PFK";
            metabolites = Dict(
                "f6p" => -1.0,
                "atp" => -1.0,
                "adp" => 1.0,
                "h" => 1.0,
                "fdp" => 1.0,
            ),
        ),
        Reaction("FBA"; metabolites = Dict("fdp" => -1.0, "dhap" => 1.0, "g3p" => 1.0)),
        Reaction("TPI"; metabolites = Dict("dhap" => -1.0, "g3p" => 1.0)),
        Reaction(
            "GAPD";
            metabolites = Dict(
                "g3p" => -1.0,
                "nad" => -1.0,
                "pi" => -1.0,
                "h" => 1.0,
                "nadh" => 1.0,
                "13dpg" => 1.0,
            ),
        ),
        Reaction(
            "PGK";
            metabolites = Dict("13dpg" => -1.0, "adp" => -1.0, "atp" => 1.0, "3pg" => 1.0),
        ),
        Reaction("PGM"; metabolites = Dict("3pg" => -1.0, "2pg" => 1)),
        Reaction("ENO"; metabolites = Dict("2pg" => -1.0, "h2o" => 1.0, "pep" => 1)),
        Reaction(
            "PYK";
            metabolites = Dict(
                "pep" => -1.0,
                "adp" => -1.0,
                "h" => -1.0,
                "atp" => 1.0,
                "pyr" => 1.0,
            ),
        ),
        Reaction(
            "LDH";
            metabolites = Dict(
                "pyr" => -1.0,
                "nadh" => -1.0,
                "h" => -1.0,
                "nad" => 1.0,
                "lac__D" => 1.0,
            ),
        ),
    ]

    add_metabolites!(model, mets)
    add_reactions!(model, rxns)

    thermodynamic_data = Dict(
        "TPI" => 5.57535,
        "PGK" => -19.32,
        "PFK" => -14.5988,
        "ENO" => -3.81089,
        "PYK" => -27.5833,
        "LDH" => -23.6803,
        "FBA" => 22.3932,
        "PGI" => 2.6617,
        "GAPD" => 4.60271,
        "PGM" => -4.52041,
        "HEX" => -17.90,
    )

    df, dgs, concens = max_min_driving_force(
        model,
        thermodynamic_data,
        Tulip.Optimizer;
        proton_id = "h",
        water_id = "h2o",
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 500)],
        concentration_ratios = [("atp", "adp", 10.0), ("nadh", "nad", 0.1)],
        constant_concentrations = [("pi", 10e-3)],
        concentration_lb = 1e-6,
        concentration_ub = 10e-3,
    )

    @test df ≈ 2.1218911811519274
end
