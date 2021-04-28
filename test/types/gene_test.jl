@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotations = Dict("sboterm" => ["sbo"], "ncbigene" => ["ads", "asds"])

    @test all(
        contains.(
            sprint(show, MIME("text/plain"), g),
            ["gene1", "gene_name", "blah", "asds"],
        ),
    )

    g2 = Gene("gene2")
    g3 = Gene("g3")
    g3.annotations = Dict("ncbigene" => ["sbo"], "ncbigene" => ["ads", "asds"])

    gd = OrderedDict(g.id => g for g in [g, g2])
    id = check_duplicate_annotations(g3, gd)
    @test id == "gene1"

    g4 = Gene("g4")
    g4.annotations = Dict("ncbigene" => ["sbo"], "ncbigene" => ["ads22", "asd22s"])
    id = check_duplicate_annotations(g4, gd)
    @test isnothing(id)
end
