
@testset "SBML import and conversion" begin
    sbmlm = load_sbml_model(model_paths["ecoli_core_model.xml"])
    m = convert(CoreModel, sbmlm)

    @test size(stoichiometry(sbmlm)) == (92, 95)
    @test size(stoichiometry(m)) == (n_metabolites(sbmlm), n_reactions(sbmlm))
    @test length(m.S.nzval) == 380
    @test length.(bounds(sbmlm)) == (95, 95)
    @test length.(bounds(m)) == (95, 95)
    @test all([length(m.xl), length(m.xu), length(m.c)] .== 95)

    @test metabolites(m)[1:3] == ["M_succoa_c", "M_ac_c", "M_fru_b"]
    @test reactions(m)[1:3] == ["R_EX_fum_e", "R_ACONTb", "R_GLNS"]

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
