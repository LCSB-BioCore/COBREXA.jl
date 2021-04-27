
sbmlfile = joinpath("data", "ecoli_core.xml")
download_data_file(
    "http://systemsbiology.ucsd.edu/sites/systemsbiology.ucsd.edu/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.xml",
    sbmlfile,
    "78692f8509fb36534f4f9b6ade23b23552044f3ecd8b48d84d484636922ae907",
)

@testset "SBML import and conversion" begin
    sbmlm = load_sbml_model(sbmlfile)
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
