using JET

function test_jet()
    @testset "JET checks" begin
        rep = report_package("BPGates";
#            ignored_modules=(
#                AnyFrameModule(Graphs.LinAlg),
#                AnyFrameModule(Graphs.SimpleGraphs),
#                AnyFrameModule(ArrayInterface),
#                AnyFrameModule(Static),
#            )
        )
        @show rep
        @test length(JET.get_reports(rep)) == 0
        #= TODO These false positives appear. Figure out how to filter them out.
        =#
    end
end

test_jet()
