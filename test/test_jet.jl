@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using BPGates

    rep = JET.report_package(BPGates; target_modules=(BPGates,))
    println(rep)
    @test length(JET.get_reports(rep)) == 0
end
