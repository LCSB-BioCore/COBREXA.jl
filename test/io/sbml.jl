
@testset "SBML import and conversion" begin
    sbmlm = load_sbml_model(model_paths["ecoli_core_model.xml"])
    m = convert(MatrixModel, sbmlm)

    @test size(stoichiometry(sbmlm)) == (92, 95)
    @test size(stoichiometry(m)) == (metabolite_count(sbmlm), variable_count(sbmlm))
    @test length(m.S.nzval) == 380
    @test length.(variable_bounds(sbmlm)) == (95, 95)
    @test length.(variable_bounds(m)) == (95, 95)
    @test all([length(m.xl), length(m.xu), length(m.c)] .== 95)

    @test metabolite_ids(m)[1:3] == ["M_13dpg_c", "M_2pg_c", "M_3pg_c"]
    @test reaction_ids(m)[1:3] == ["R_ACALD", "R_ACALDt", "R_ACKr"]

    cm = convert(MatrixModelWithCoupling, sbmlm)
    @test n_coupling_constraints(cm) == 0
end

@testset "Save SBML model" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.xml"])
    testpath = tmpfile("modeltest.xml")
    save_model(convert(SBMLModel, model), testpath)
    wrote = convert(MatrixModel, load_sbml_model(testpath))
    @test isequal(model, wrote)
end

@testset "Import yeast-GEM (sbml)" begin
    m = load_model(ObjectModel, model_paths["yeast-GEM.xml"])
    @test metabolite_count(m) == 2744
    @test reaction_count(m) == 4063
    @test n_genes(m) == 1160
end
