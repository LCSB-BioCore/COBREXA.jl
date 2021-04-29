@testset "Metabolite" begin
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ads", "asds"])

    @test all(
        contains.(
            sprint(show, MIME("text/plain"), m1),
            ["met1", "metabolite 1", "C6H12O6N", "blah", "asds"],
        ),
    )

    m2 = Metabolite("met2")

    m2.formula = "C6H12O6N"

    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ad2s", "asds"])

    ats = get_atoms(m1)
    @test ats["C"] == 6 && ats["N"] == 1

    m4 = Metabolite("met4")
    m4.formula = "X"
    m4.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["adxxx2s", "asdxxxs"])

    md = OrderedDict(m.id => m for m in [m1, m2, m3])
    id = check_duplicate_annotations(m4, md)
    @test isnothing(id)
end
