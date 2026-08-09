@testitem "Aqua analysis" tags=[:aqua] begin

using Aqua, BPGates

Aqua.test_all(BPGates)

end
