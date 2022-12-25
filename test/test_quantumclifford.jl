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

    state = rand(BellState,num_bell)
    circuit1 = [rand(BellGate,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    circuit2 = [rand(BellPauliPermutation,randperm(num_bell)[1:1]...) for _ in 1:num_gates]
    circuit3 = [rand(BellSinglePermutation,randperm(num_bell)[1:1]...) for _ in 1:num_gates]
    circuit4 = [rand(BellDoublePermutation,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    # TODO implement CNOTPerm conversions
    circuit5 = []#[rand(CNOTPerm,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    circuit6 = [CNOTPerm(1,1,randperm(num_bell)[1:2]...) for _ in 1:num_gates]
    circuit = [circuit1...;circuit2...;circuit3...;circuit4...;circuit5...;circuit6...;]
    endstate, status = mctrajectory!(copy(state), circuit)

    stabstate = MixedDestabilizer(copy(state))
    stabcircuit = [toQCcircuit.(circuit)...;]
    endstabstate, stabstatus = mctrajectory!(copy(stabstate), stabcircuit)

    @test status == stabstatus
    @test canonicalize!(copy(stabilizerview(endstabstate))) == canonicalize!(Stabilizer(endstate))

    # Test measurements separately because resets to Bell states can get shortcircuited
    meas = rand(BellMeasure,1)
    meas_qc = [toQCcircuit(meas)...]
    mstate, mstatus = mctrajectory!(endstate, [meas])
    mstate_qc, mstatus_qc = mctrajectory!(endstabstate, meas_qc)
    mstate_qc = apply!(mstate_qc, Reset(S"XX ZZ",[1,2]))
    @test mstatus == mstatus_qc
    @test canonicalize!(copy(stabilizerview(mstate_qc))) == canonicalize!(Stabilizer(mstate))
end

state = BellState(2)
for _ in 1:10
    new_state = apply!(copy(state), rand(CNOTPerm,1,2))
    @test state == new_state
end
