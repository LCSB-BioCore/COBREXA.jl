@testset "Gene" begin
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes"=>["blah", "blah"])
    g.annotation = Dict("sboterm" => "sbo", "ncbigene" => ["ads", "asds"])

    @test repr("text/plain", g) == "Gene ID: gene1\nGene name: gene_name\n"

    g2 = Gene("gene2")
   
    genes = [g, g2]
    @test repr("text/plain", genes) == "Gene set of length: 2\n"

    gene_list = [[g], [g2]]
    @test repr("text/plain", gene_list) == "(gene1) or (gene2)\n"

    @test genes[g] == 1
    
    gg = findfirst(genes, g2.id)
    @test gg.id == g2.id

    g3 = Gene("g3")
    g3.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads", "asds"])
    
    dup, ind = check_duplicate_annotations(genes, g3)
    @test dup && ind == 1
end
