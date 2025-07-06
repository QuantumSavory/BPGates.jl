@testitem "Quantikz circuit plotting" begin

using BPGates
using Quantikz
using Test

@testset "Quantikz" begin
    quantikz_circ = Quantikz.QuantikzOp.([BellSinglePermutation(2,2), BellDoublePermutation(1,2,3), BellPauliPermutation(4,1), BellMeasure(3,1), BellGate(1,1,1,1,1,1,2), CNOTPerm(1,2,3,4), GoodSingleQubitPerm(1,2), PauliNoiseOp(1,0.0,0.0,0.0), PauliNoiseBellGate(CNOTPerm(1,2,3,4),0.0,0.0,0.0), NoisyBellMeasure(BellMeasure(1,2),0.0), NoisyBellMeasureNoisyReset(BellMeasure(1,2),0.0,0.0,0.0,0.0),T1NoiseOp(1,0.1),T2NoiseOp(2,0.2)])
    img = Quantikz.circuit2image(quantikz_circ)
end

end
