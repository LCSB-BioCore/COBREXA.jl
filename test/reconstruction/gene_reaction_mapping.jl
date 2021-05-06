@testset "mapping 1 gene:1 reaction" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Reaction("v1", metabolites = Dict{String,Int}(), grr = [["g1"]]))
    @test m.genes["g1"].associated_reactions == Set(("v1",))

    rm!(Reaction, m, "v1")
    @test m.genes["g1"].associated_reactions == Set{String}()
end

@testset "mapping 1 gene:N reactions" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites = Dict{String,Int}(), grr = [["g1"]]))
    add!(m, Reaction("v2", metabolites = Dict{String,Int}(), grr = [["g1"]]))
    @test m.genes["g1"].associated_reactions == Set(("v1", "v2"))

    rm!(Reaction, m, "v1")
    @test m.genes["g1"].associated_reactions == Set(("v2",))
end

@testset "mapping 1 genes:1 reaction" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites = Dict{String,Int}(), grr = [["g1"], ["g2"]]))
    @test m.genes["g1"].associated_reactions == Set(("v1",))
    @test m.genes["g2"].associated_reactions == Set(("v1",))

    rm!(Reaction, m, "v1")
    @test m.genes["g1"].associated_reactions == Set{String}()
    @test m.genes["g1"].associated_reactions == Set{String}()
end

@testset "mapping N genes:N reactions" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites = Dict{String,Int}(), grr = [["g1"], ["g2"]]))
    add!(m, Reaction("v2", metabolites = Dict{String,Int}(), grr = [["g1"], ["g2"]]))
    @test m.genes["g1"].associated_reactions == Set(("v1", "v2"))
    @test m.genes["g2"].associated_reactions == Set(("v1", "v2"))

    rm!(Reaction, m, "v1")
    @test m.genes["g1"].associated_reactions == Set(("v2",))
    @test m.genes["g2"].associated_reactions == Set(("v2",))
end
