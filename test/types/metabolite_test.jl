@testset "Metabolite" begin
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotation = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), m1) ==
          "\e[34mMetabolite ID: \e[35mmet1\n\e[34mName: \e[35mmetabolite 1\n\e[34mFormula: \e[35mC6H12O6N\n\e[34mCharge: \e[35m1\n\e[34mCompartment: \e[35mc\n\e[34mNotes: \n\e[35m\tnotes: blah, blah\n\e[34mAnnotation: \n\e[35m\tkegg.compound: ads, asds\n\e[35m\tsboterm: sbo\n\e[34mFields: \e[35mid, name, formula, charge, compartment, notes, annotation\n"

    m2 = Metabolite("met2")
    m2.formula = "C6H12O6N"

    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotation = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ad2s", "asds"])

    mets = [m1, m2, m3]

    @test sprint(show, MIME("text/plain"), mets) ==
          "\e[34mMetabolite vector of length: : \e[35m3\n\e[34mEach metabolite has fields: \e[35mid, name, formula, charge, compartment, notes, annotation\n"

    md = OrderedDict(m.id => m for m in mets)
    id = check_duplicate_annotations(m3, md)
    @test id == "met3"

    ats = get_atoms(m1)
    @test ats["C"] == 6 && ats["N"] == 1

    m4 = Metabolite("met4")
    m4.formula = "X"
    m4.annotation = Dict("sboterm" => ["sbo"], "kegg.compound" => ["adxxx2s", "asdxxxs"])

    id = check_duplicate_annotations(m4, md)
    @test isnothing(id)
end
