@testset "Gene: construction, printing, utils" begin
    g = Gene()

    # test defaults
    @test isnothing(g.name)
    @test isempty(g.notes)
    @test isempty(g.annotations)
    
    # Now assign
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotations = Dict("sboterm" => ["sbo"], "ncbigene" => ["ads", "asds"])

    # Test pretty printing
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), g),
            ["gene1", "gene_name", "blah", "asds"],
        ),
    )

    # Test duplicate annotation finder
    g2 = Gene("gene2")
    g2.annotations = Dict("sboterm" => ["sbo2"], "ncbigene" => ["fff", "ggg"])
    g3 = Gene("g3")
    g3.annotations = Dict("sboterm" => ["sbo3"], "ncbigene" => ["ads"])
    g4 = Gene("g4")
    g4.annotations = Dict("sboterm" => ["sbo4"], "ncbigene" => ["ads22", "asd22s"])
    gdict = OrderedDict(g.id => g for g in [g, g2]) # this is how genes are stored in StandardModel

    id = check_duplicate_annotations(g3, gdict)
    @test id == "gene1"

    id = check_duplicate_annotations(g4, gdict)
    @test isnothing(id)
end
