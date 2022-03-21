@testset "Construction" begin
    cm = CommunityModel()
    @test isa(cm, CommunityModel)
    cm = CommunityModel(test_toyModel())
    @test isa(cm, CommunityModel)
end

@testset "Basic getters" begin
    cm = CommunityModel(test_toyModel())
    @test reactions(cm) == reactions(test_toyModel())
    @test metabolites(cm) == metabolites(test_toyModel())
    @test stoichiometry(cm) == stoichiometry(test_toyModel())
    cm = CommunityModel(test_LP())
    @test bounds(cm) == bounds(test_LP())
    @test objective(cm) == objective(test_LP())
end
