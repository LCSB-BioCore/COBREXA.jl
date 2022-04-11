@testset "SMOMENT" begin
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

    remove_slow_isozymes!(model, rid_isozymes)

    smm = SMomentModel(model; rid_isozymes, enzyme_capacity = total_protein_mass)

    change_bounds(
        smm,
        ["EX_glc__D_e", "GLCpts"];
        lower = [-1000.0, -1.0],
        upper = [nothing, 12.0],
    )

    rxn_fluxes = flux_balance_analysis_dict(
        smm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.8907273630431708,
        atol = TEST_TOLERANCE,
    )
end
