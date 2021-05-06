COBREXA.@_make_logging_tag TEST "testing stuff"

log_TEST()
@testset "Logging on" begin
    @test _TEST_log_enabled
    @test_logs (:warn, "qweasdzxc") @_TEST_log @warn "qweasdzxc"
end

log_TEST(false)
@testset "Logging off" begin
    @test !_TEST_log_enabled
    @test_logs @_TEST_log @warn "all okay!"
end
