@testset "Metabolite" begin
    m1 = Metabolite(id="testmetabolite")
    m1.id = "met1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ads", "asds"])

    @test all(
        contains.(
            sprint(show, MIME("text/plain"), m1),
            ["met1", "C6H12O6N", "blah", "asds"],
        ),
    )

    m2 = Metabolite(id="met2")

    m2.formula = "C6H12O6N"

    m3 = Metabolite(id="met3")
    m3.formula = "X"
    m3.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ad2s", "asds"])

    m4 = Metabolite(id="met4")
    m4.formula = "X"
    m4.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["adxxx2s", "asdxxxs"])

    md = OrderedDict(m.id => m for m in [m1, m2, m3])
    @test issetequal(["met1", "met3"], ambiguously_identified_items(annotation_index(md)))
end
