@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotation = Dict("sboterm" => ["sbo"], "ncbigene" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), g) == "Gene.id: gene1\nGene.name: gene_name\nGene.notes: \n\tnotes: blah, blah\nGene.annotation: \n\tncbigene: ads, asds\n\tsboterm: sbo\n"
    
    g2 = Gene("gene2")

    genes = [g, g2]

    gene_list = [[g], [g2]]
    @test sprint(show, MIME("text/plain"), gene_list) == "Gene reaction rule: (gene1) or (gene2)\n"

    @test genes[g] == 1

    gg = findfirst(genes, g2.id)
    @test gg.id == g2.id

    g3 = Gene("g3")
    g3.annotation = Dict("ncbigene" => ["sbo"], "ncbigene" => ["ads", "asds"])

    gd = OrderedDict(g.id => g for g in genes)
    id = check_duplicate_annotations(g3, gd)
    @test id == "gene1"

    g4 = Gene("g4")
    g4.annotation = Dict("ncbigene" => ["sbo"], "ncbigene" => ["ads22", "asd22s"])
    id = check_duplicate_annotations(g4, gd)
    @test isnothing(id)
end
