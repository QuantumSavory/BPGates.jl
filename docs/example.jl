using Revise
using BenchmarkTools
using BPGates, QuantumClifford

num_bell = 4

##################
## apply gates
##################

gate = BellGateQC(1,14,2,1,1,2)
# range of inputs: pauli: 1-4, double: 1-20, single1: 1-6, single2: 1-6, index1: 1-num_bell, index2: 1-num_bell
gate1 = BellPauliPermutation(1,(1,2))
gate2 = BellDoublePermutation(14,(1,2))
gate3 = BellSinglePermutation(2,1)

state = rand_state(num_bell)
state2 = copy(state) # copy the state

# @benchmark apply!(state, gate)
# @benchmark apply_as_qc!(state2, gate)
apply!(state, gate)
apply_as_qc!(state2, gate)

@assert state.phases == state2.phases # resultant state after applying gates are the same

res = bellmeasure!(state, BellMeasure(1,1)) # measure index: 1-3, bell index: 1-num_bell
res2 = apply_as_qc!(state2, BellMeasure(1,1))

@assert res[end] == res2[end] # measurement results are the same

##################
## apply circuits
##################
num_bell = 8

state = rand_state(num_bell)
state2 = copy(state) # copy the state

circuit = [BellGateQC(1,14,2,1,1,2), BellGateQC(1,13,2,2,1,2)]

# @benchmark QuantumClifford.Experimental.NoisyCircuits.mctrajectory!(state, circuit)
# @benchmark apply_as_qc!(state2, circuit)
QuantumClifford.Experimental.NoisyCircuits.mctrajectory!(state, circuit)
apply_as_qc!(state2, circuit)

@assert state.phases == state2.phases # resultant state after applying gates are the same

res = bellmeasure!(state, BellMeasure(1,1)) # measure index: 1-3, bell index: 1-num_bell
res2 = apply_as_qc!(state2, BellMeasure(1,1))

@assert res[end] == res2[end] # measurement results are the same
