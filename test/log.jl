COBREXA.Log.Internal.@make_logging_tag TEST "testing stuff"

log_TEST()
@testset "Logging on" begin
    @test TEST_log_enabled
    @test_logs (:warn, "qweasdzxc") @TEST_log @warn "qweasdzxc"
end

log_TEST(false)
@testset "Logging off" begin
    @test !TEST_log_enabled
    @test_logs @TEST_log @warn "all okay!"
end
