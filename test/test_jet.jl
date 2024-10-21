@testitem "JET analysis" tags=[:jet] begin

using JET
using Test
using BPGates

rep = report_package("BPGates";
    ignored_modules=(
        LastFrameModule(Base),
    )
)
@show rep
@test length(JET.get_reports(rep)) == 0

end
