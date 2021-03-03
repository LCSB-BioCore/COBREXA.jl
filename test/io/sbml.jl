
sbmlfile = joinpath("data", "ecoli_core.xml")
if !isfile(sbmlfile)
    download(
        "http://systemsbiology.ucsd.edu/sites/systemsbiology.ucsd.edu/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.xml",
        sbmlfile,
    )
end

cksum = bytes2hex(sha256(open(sbmlfile)))
if cksum != "78692f8509fb36534f4f9b6ade23b23552044f3ecd8b48d84d484636922ae907"
    @warn "The downloaded E Coli core model seems to be different from the expected one. Tests may fail." cksum
end

@testset "SBML import" begin
    m = loadSBMLModel(sbmlfile)

    @test size(m.C) == (0, 95)
    @test size(m.S) == (92, 95)
    @test length(m.S.nzval) == 380
    @test length(m.b) == 92
    @test all([length(m.xl), length(m.xu), length(m.c)] .== 95)

    @test m.mets[1:3] == ["M_succoa_c", "M_ac_c", "M_fru_b"]
    @test m.rxns[1:3] == ["R_EX_fum_e", "R_ACONTb", "R_GLNS"]
end
