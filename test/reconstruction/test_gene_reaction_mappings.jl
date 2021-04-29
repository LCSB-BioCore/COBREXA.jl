@testset "single gene - single reaction" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Reaction("v1", metabolites=Dict{String, Int}(), grr=[["g1"]]))
    @test m.genes["g1"].reactions == Set(("v1", ))
    
    rm!(Reaction, m, "v1")
    @test m.genes["g1"].reactions == Set{String}()
end

@testset "single gene - multiple reactions" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites=Dict{String, Int}(), grr=[["g1"]]))
    add!(m, Reaction("v2", metabolites=Dict{String, Int}(), grr=[["g1"]]))
    @test m.genes["g1"].reactions == Set(("v1", "v2"))
    
    rm!(Reaction, m, "v1")
    @test m.genes["g1"].reactions == Set(("v2", ))
end

@testset "multiple genes - single reactions" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites=Dict{String, Int}(), grr=[["g1"], ["g2"]]))
    @test m.genes["g1"].reactions == Set(("v1",))
    @test m.genes["g2"].reactions == Set(("v1",))
    
    rm!(Reaction, m, "v1")
    @test m.genes["g1"].reactions == Set{String}()
    @test m.genes["g1"].reactions == Set{String}()
end

@testset "multiple genes - multiple reactions" begin
    m = StandardModel()
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites=Dict{String, Int}(), grr=[["g1"], ["g2"]]))
    add!(m, Reaction("v2", metabolites=Dict{String, Int}(), grr=[["g1"], ["g2"]]))
    @test m.genes["g1"].reactions == Set(("v1", "v2"))
    @test m.genes["g2"].reactions == Set(("v1", "v2"))
    
    rm!(Reaction, m, "v1")
    @test m.genes["g1"].reactions == Set(("v2", ))
    @test m.genes["g2"].reactions == Set(("v2", ))
end
