using JET

function test_jet()
    @testset "JET checks" begin
        rep = report_package("BPGates";
#            ignored_modules=(
#            )
        )
        @show rep
        @test length(JET.get_reports(rep)) == 0
    end
end

test_jet()
