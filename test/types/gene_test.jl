@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotation = Dict("sboterm" => "sbo", "ncbigene" => ["ads", "asds"])

    @test sprint(show, MIME("text/plain"), g) ==
          "\e[34mGene.id: \e[35mgene1\n\e[34mGene.name: \e[35mgene_name\n\e[34mGene.notes: \n\e[35m\tnotes: blah, blah\n\e[34mGene.annotation: \n\e[35m\tncbigene: ads, asds\n\e[35m\tsboterm: sbo\n"

    g2 = Gene("gene2")

    genes = [g, g2]

    gene_list = [[g], [g2]]
    @test sprint(show, MIME("text/plain"), gene_list) ==
          "\e[34mGene reaction rule: \e[35m(gene1) or (gene2)\n"

    @test genes[g] == 1

    gg = findfirst(genes, g2.id)
    @test gg.id == g2.id

    g3 = Gene("g3")
    g3.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads", "asds"])

    dup, ind = check_duplicate_annotations(genes, g3)
    @test dup && ind == 1

    @test isnothing(findfirst(genes, "nope"))

    g4 = Gene("g4")
    g4.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads22", "asd22s"])
    dup, ind = check_duplicate_annotations(genes, g4)
    @test !dup && ind == -1
end
