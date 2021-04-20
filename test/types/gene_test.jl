@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotation = Dict("sboterm" => ["sbo"], "ncbigene" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), g) == "\e[34mGene ID: \e[35mgene1\n\e[34mName: \e[35mgene_name\n\e[34mNotes: \n\e[35m\tnotes: blah, blah\n\e[34mAnnotation: \n\e[35m\tncbigene: ads, asds\n\e[35m\tsboterm: sbo\n\e[34mFields: \e[35mid, name, notes, annotation\n"

    g2 = Gene("gene2")

    genes = [g, g2]
    @test sprint(show, MIME("text/plain"), genes) == "\e[34mGene vector of length: \e[35m2\n\e[34mEach gene has fields: \e[35mid, name, notes, annotation\n"

    gene_list = [[g], [g2]]
    @test sprint(show, MIME("text/plain"), gene_list) == "\e[34mGene reaction rule: \e[35m(gene1) or (gene2)\n"

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
