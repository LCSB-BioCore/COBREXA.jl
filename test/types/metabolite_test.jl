@testset "Metabolite" begin
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), m1) ==
          "Metabolite ID: met1\nMetabolite name: metabolite 1\nFormula: C6H12O6N\nCharge: 1\n"

    m2 = Metabolite("met2")
    m2.formula = "C6H12O6N"

    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ad2s", "asds"])

    mets = [m1, m2, m3]

    @test sprint(show, MIME("text/plain"), mets) == "Metabolite set of length: 3\n"

    mets_d = OrderedDict(zip([m.id for m in mets], mets))
    dup, ind = check_duplicate_annotations(m3, mets_d)
    @test dup && ind == "met3"

    ats = get_atoms(mets[1])
    @test ats["C"] == 6 && ats["N"] == 1

    m4 = Metabolite("met4")
    m4.formula = "X"
    m4.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["adxxx2s", "asdxxxs"])

    dup, ind = check_duplicate_annotations(m4, mets_d)
    @test !dup && ind == ""
end
