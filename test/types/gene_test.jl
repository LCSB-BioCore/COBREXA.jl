@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotation = Dict("sboterm" => "sbo", "ncbigene" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), g) == "Gene ID: gene1\nGene name: gene_name\n"

    g2 = Gene("gene2")

    genelist = [g, g2]
    @test sprint(show, MIME("text/plain"), genelist) == "Gene set of length: 2\n"

    g3 = Gene("g3")
    g3.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads", "asds"])

    genes_d = OrderedDict(zip([g.id for g in genelist], genelist))
    dup, ind = check_duplicate_annotations(g3, genes_d)
    @test dup && ind == "gene1"

    g4 = Gene("g4")
    g4.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads22", "asd22s"])
    dup, ind = check_duplicate_annotations(g4, genes_d)
    @test !dup && ind == ""
end
