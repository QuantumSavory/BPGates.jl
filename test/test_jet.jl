using JET
using ArrayInterface
using Static

function test_jet()
    @testset "JET checks" begin
        rep = report_package("BPGates";
            ignored_modules=(
                AnyFrameModule(ArrayInterface),
                AnyFrameModule(Static),
            )
        )
        @show rep
        @test length(JET.get_reports(rep)) == 0
    end
end

test_jet()
