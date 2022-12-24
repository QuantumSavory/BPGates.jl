using BPGates, QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using Random
using BPGates: toQCcircuit
using Test

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

# TODO finish and move to the actual library
function BPGates.toQCcircuit(g::CNOTPerm)
    return [
        # MISSING STEPS
        sCNOT(g.idx1*2-1, g.idx2*2-1),
        sCNOT(g.idx1*2, g.idx2*2)
    ]
end

for num_bell in test_sizes[2:end]
    num_gates = 10

    state = BellState(num_bell)
    circuit1 = [rand(BellGate,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    circuit2 = [rand(BellPauliPermutation,randperm(num_bell)[1:1]...) for _ in 1:num_gates]
    circuit3 = [rand(BellSinglePermutation,randperm(num_bell)[1:1]...) for _ in 1:num_gates]
    circuit4 = [rand(BellDoublePermutation,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    # TODO implement CNOTPerm conversions
    circuit5 = []#[rand(CNOTPerm,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    circuit6 = [CNOTPerm(1,1,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    # TODO implement BellMeasure conversions
    circuit7 = []#[rand(BellMeasure,i) for i in 1:num_bell÷3]
    circuit = [circuit1...;circuit2...;circuit3...;circuit4...;circuit5...;circuit6...;circuit7...;]
    endstate, status = mctrajectory!(copy(state), circuit)

    stabstate = Stabilizer(copy(state))
    stabcircuit = [toQCcircuit.(circuit)...;]
    endstabstate, stabstatus = mctrajectory!(copy(stabstate), stabcircuit)

    @test canonicalize!(endstabstate) == canonicalize!(Stabilizer(endstate))
    @test status == stabstatus
end

state = BellState(2)
for _ in 1:10
    new_state = apply!(copy(state), rand(CNOTPerm,1,2))
    @test state == new_state
end
