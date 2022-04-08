@testset "GECKO" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])
    total_protein_mass = 100 # mg/gdW

    #: construct isozymes from model
    rid_isozymes = Dict{String,Vector{Isozyme}}()
    for (rid, kcats) in ecoli_core_reaction_kcats
        grrs = reaction_gene_association(model, rid)
        rid_isozymes[rid] = [
            Isozyme(
                Dict(grrs[i] .=> ecoli_core_protein_stoichiometry[rid][i]),
                (kcats[i][1], kcats[i][2]),
            ) for i = 1:length(grrs)
        ]
    end

    #: add molar mass to genes in model
    for (gid, g) in model.genes
        model.genes[gid].molar_mass = get(ecoli_core_protein_masses, gid, nothing)
    end

    gm = GeckoModel(
        model;
        rid_isozymes,
        enzyme_capacities = [(get_genes_with_kcats(rid_isozymes), total_protein_mass)],
    )
    change_bounds(
        gm,
        ["EX_glc__D_e", "b2779", "GLCpts"];
        lbs = [-1000.0, 0.01, -1.0],
        ubs = [nothing, 0.06, 12.0],
    )

    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
        sense = COBREXA.MOI.MAX_SENSE,
    )

    rxn_fluxes = flux_dict(gm, opt_model)
    prot_concens = protein_dict(gm, opt_model)

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.812827846796761,
        atol = TEST_TOLERANCE,
    )

    prot_mass = sum(ecoli_core_protein_masses[gid] * c for (gid, c) in prot_concens)

    @test isapprox(prot_mass, total_protein_mass, atol = TEST_TOLERANCE)
end

