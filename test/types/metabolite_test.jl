@testset "Metabolite" begin
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), m1) == "\e[34mMetabolite.id: \e[35mmet1\n\e[34mMetabolite.name: \e[35mmetabolite 1\n\e[34mMetabolite.formula: \e[35mC6H12O6N\n\e[34mMetabolite.charge: \e[35m1\n\e[34mMetabolite.compartment: \e[35mc\n\e[34mMetabolite.notes: \n\e[35m\tnotes: blah, blah\n\e[34mMetabolite.annotation: \n\e[35m\tkegg.compound: ads, asds\n\e[35m\tsboterm: sbo\n"

    m2 = Metabolite("met2")

    m2.formula = "C6H12O6N"

    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ad2s", "asds"])

    mets = [m1, m2, m3]

    @test mets[m2] == 2

    mm = findfirst(mets, "met3")
    @test mm.id == m3.id

    dup, ind = check_duplicate_annotations(mets, m3)
    @test dup && ind == 3

    mms = check_same_formula([m3, m1], m2)
    @test length(mms) == 1

    ats = get_atoms(mms[1])
    @test ats["C"] == 6 && ats["N"] == 1

    @test isnothing(findfirst(mets, "nope"))

    m4 = Metabolite("met4")
    m4.formula = "X"
    m4.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["adxxx2s", "asdxxxs"])

    dup, ind = check_duplicate_annotations(mets, m4)
    @test !dup && ind == -1
end
