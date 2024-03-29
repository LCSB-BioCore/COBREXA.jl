
@testset "SBML import and conversion" begin
    sbmlm = load_sbml_model(model_paths["ecoli_core_model.xml"])
    m = convert(CoreModel, sbmlm)

    @test size(stoichiometry(sbmlm)) == (92, 95)
    @test size(stoichiometry(m)) == (n_metabolites(sbmlm), n_reactions(sbmlm))
    @test length(m.S.nzval) == 380
    @test length.(bounds(sbmlm)) == (95, 95)
    @test length.(bounds(m)) == (95, 95)
    @test all([length(m.xl), length(m.xu), length(m.c)] .== 95)

    @test metabolites(m)[1:3] == ["M_13dpg_c", "M_2pg_c", "M_3pg_c"]
    @test reactions(m)[1:3] == ["R_ACALD", "R_ACALDt", "R_ACKr"]

    cm = convert(CoreModelCoupled, sbmlm)
    @test n_coupling_constraints(cm) == 0
end

@testset "Save SBML model" begin
    model = load_model(CoreModel, model_paths["e_coli_core.xml"])
    testpath = tmpfile("modeltest.xml")
    save_model(convert(SBMLModel, model), testpath)
    wrote = convert(CoreModel, load_sbml_model(testpath))
    @test isequal(model, wrote)
end

@testset "Import yeast-GEM (sbml)" begin
    m = load_model(StandardModel, model_paths["yeast-GEM.xml"])
    @test n_metabolites(m) == 2744
    @test n_reactions(m) == 4063
    @test n_genes(m) == 1160
end
