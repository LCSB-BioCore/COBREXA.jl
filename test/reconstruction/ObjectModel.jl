@testset "Model manipulation" begin
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    mets = [m1, m2, m3, m4]
    m5 = Metabolite("m5")
    m6 = Metabolite("m6")
    m7 = Metabolite("m7")
    mtest = Metabolite("mtest")

    g1 = Gene("g1")
    g2 = Gene("g2")
    g3 = Gene("g3")
    g4 = Gene("g4")
    gene_vec = [g1, g2, g3, g4]
    g5 = Gene("g5")
    g6 = Gene("g6")
    g7 = Gene("g7")
    gtest = Gene("gtest")

    r1 = ReactionForward("r1", Dict(m1.id => -1.0, m2.id => 1.0))
    r2 = ReactionBidirectional("r2", Dict(m2.id => -2.0, m3.id => 1.0))
    r2.gene_associations = [Isozyme(x) for x in [["g2"], ["g1", "g3"]]]
    r3 = ReactionBackward("r3", Dict(m1.id => -1.0, m4.id => 2.0))
    r4 = ReactionBackward("r4", Dict(m1.id => -5.0, m4.id => 2.0))
    r5 = ReactionBackward("r5", Dict(m1.id => -11.0, m4.id => 2.0, m3.id => 2.0))
    rtest = ReactionForward("rtest", Dict(m1.id => -1.0, m2.id => 1.0))

    rxns = [r1, r2]

    model = ObjectModel()
    model.reactions = OrderedDict(r.id => r for r in rxns)
    model.metabolites = OrderedDict(m.id => m for m in mets)
    model.genes = OrderedDict(g.id => g for g in gene_vec)

    # change bound tests - in place
    change_bound!(model, "r2"; lower_bound = -10, upper_bound = 10)
    @test model.reactions["r2"].lower_bound == -10
    @test model.reactions["r2"].upper_bound == 10
    change_bound!(model, "r2"; lower_bound = -100)
    @test model.reactions["r2"].lower_bound == -100
    @test model.reactions["r2"].upper_bound == 10
    change_bound!(model, "r2"; upper_bound = 111)
    @test model.reactions["r2"].lower_bound == -100
    @test model.reactions["r2"].upper_bound == 111

    change_bounds!(
        model,
        ["r1", "r2"];
        lower_bounds = [-110, -220],
        upper_bounds = [110.0, 220.0],
    )
    @test model.reactions["r1"].lower_bound == -110
    @test model.reactions["r1"].upper_bound == 110
    @test model.reactions["r2"].lower_bound == -220
    @test model.reactions["r2"].upper_bound == 220

    # change bound - new model
    new_model = change_bound(model, "r2"; lower_bound = -10, upper_bound = 10)
    @test new_model.reactions["r2"].lower_bound == -10
    @test new_model.reactions["r2"].upper_bound == 10

    new_model = change_bound(model, "r2"; lower_bound = -10)
    @test new_model.reactions["r2"].lower_bound == -10
    @test new_model.reactions["r2"].upper_bound == 220

    new_model = change_bounds(
        model,
        ["r1", "r2"];
        lower_bounds = [-10, -20],
        upper_bounds = [10.0, 20.0],
    )
    @test new_model.reactions["r1"].lower_bound == -10
    @test new_model.reactions["r1"].upper_bound == 10
    @test new_model.reactions["r2"].lower_bound == -20
    @test new_model.reactions["r2"].upper_bound == 20

    new_model = change_bounds(
        model,
        ["r1", "r2"];
        lower_bounds = [-10, nothing],
        upper_bounds = [nothing, 20.0],
    )
    @test new_model.reactions["r1"].lower_bound == -10
    @test new_model.reactions["r1"].upper_bound == 110
    @test new_model.reactions["r2"].lower_bound == -220
    @test new_model.reactions["r2"].upper_bound == 20

    ### objective
    change_objective!(model, "r2")
    @test model.objective["r2"] == 1.0

    new_model = change_objective(model, "r1"; weight = 2.0)
    @test new_model.objective["r1"] == 2.0

    ### reactions
    add_reactions!(model, [r3, r4])
    @test length(model.reactions) == 4

    add_reaction!(model, r5)
    @test length(model.reactions) == 5

    remove_reactions!(model, ["r5", "r4"])
    @test length(model.reactions) == 3

    remove_reaction!(model, "r1")
    @test length(model.reactions) == 2

    new_model = model |> with_added_reaction(rtest)
    @test haskey(new_model.reactions, "rtest")
    @test !haskey(model.reactions, "rtest")

    new_model2 = new_model |> with_removed_reaction("rtest")
    @test !haskey(new_model2.reactions, "rtest")
    @test haskey(new_model.reactions, "rtest")

    @test_throws DomainError add_reaction!(model, r3)
    @test_throws DomainError remove_reaction!(model, "abc")

    ### metabolites
    add_metabolites!(model, [m5, m6])
    @test length(model.metabolites) == 6

    add_metabolite!(model, m7)
    @test length(model.metabolites) == 7

    remove_metabolites!(model, ["m5", "m4"])
    @test length(model.metabolites) == 5

    remove_metabolite!(model, "m1")
    @test length(model.metabolites) == 4

    new_model = model |> with_added_metabolite(mtest)
    @test haskey(new_model.metabolites, "mtest")
    @test !haskey(model.metabolites, "mtest")

    new_model2 = new_model |> with_removed_metabolite("mtest")
    @test !haskey(new_model2.metabolites, "mtest")
    @test haskey(new_model.metabolites, "mtest")

    @test_throws DomainError add_metabolite!(model, m2)
    @test_throws DomainError remove_metabolite!(model, "abc")

    ### genes
    add_genes!(model, [g5, g6])
    @test length(model.genes) == 6

    add_gene!(model, g7)
    @test length(model.genes) == 7

    remove_genes!(model, ["g5", "g4"])
    @test length(model.genes) == 5

    remove_gene!(model, "g1")
    @test length(model.genes) == 4

    new_model = model |> with_added_gene(gtest)
    @test haskey(new_model.genes, "gtest")
    @test !haskey(model.genes, "gtest")

    new_model2 = new_model |> with_removed_gene("gtest")
    @test !haskey(new_model2.genes, "gtest")
    @test haskey(new_model.genes, "gtest")

    @test_throws DomainError add_gene!(model, g7)
    @test_throws DomainError remove_gene!(model, "abc")

    # change gene
    change_gene_product_bound!(model, "g3"; lower_bound = -10, upper_bound = 10)
    @test model.genes["g3"].product_lower_bound == -10.0
    @test model.genes["g3"].product_upper_bound == 10.0

    new_model = change_gene_product_bound(model, "g2"; lower_bound = -10, upper_bound = 10)
    @test new_model.genes["g2"].product_lower_bound == -10.0
    @test new_model.genes["g2"].product_upper_bound == 10.0
    @test model.genes["g2"].product_lower_bound == 0.0

    # isozymes
    isos = [Isozyme(["g1"])]
    new_model = model |> with_removed_isozymes("r2")
    @test isnothing(new_model.reactions["r2"].gene_associations)
    @test !isnothing(model.reactions["r2"].gene_associations)

    new_model2 = new_model |> with_added_isozymes("r2", isos)
    @test !isnothing(new_model2.reactions["r2"].gene_associations)
    @test isnothing(new_model.reactions["r2"].gene_associations)
end
