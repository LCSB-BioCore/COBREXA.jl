@testset "Metabolite" begin
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotation = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), m1) == "Metabolite ID: met1\nName: metabolite 1\nFormula: C6H12O6N\nCharge: 1\nCompartment: c\nNotes: \n\tnotes: blah, blah\nAnnotation: \n\tkegg.compound: ads, asds\n\tsboterm: sbo\nFields: id, name, formula, charge, compartment, notes, annotation\n"

    m2 = Metabolite("met2")

    m2.formula = "C6H12O6N"

    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotation = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ad2s", "asds"])

    mets = [m1, m2, m3]

    @test sprint(show, MIME("text/plain"), mets) == "Metabolite vector of length: : 3\nEach metabolite has fields: id, name, formula, charge, compartment, notes, annotation\n"

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
