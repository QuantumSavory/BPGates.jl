using Revise
using BenchmarkTools
using BPGates
using QuantumClifford

num_bell = 4

##################
## apply gates
##################

state = rand_state(num_bell)
state2 = copy(state) # copy the state

gate = BellGateQC(1,14,2,1,1,2)
# range of inputs: pauli: 1-6, double: 1-20, single1: 1-6, single2: 1-6, index1: 1-num_bell, index2: 1-num_bell

apply_op!(state, gate)
apply_as_qc!(state2, gate)

@assert state.phases == state2.phases # resultant state after applying gates are the same

res = apply_op!(state, BellMeasure(1,1))[end] # measure index: 1-3, bell index: 1-num_bell
res2 = apply_as_qc!(state2, BellMeasure(1,1))[end]

@assert res == res2 # measurement results are the same

##################
## apply circuits
##################

state = rand_state(num_bell)
state2 = copy(state) # copy the state

circuit = [BellGateQC(1,14,2,1,1,2), BellGateQC(1,13,2,2,1,2)]

apply_op!(state, circuit)
apply_as_qc!(state2, circuit)

@assert state.phases == state2.phases # resultant state after applying gates are the same

res = apply_op!(state, BellMeasure(1,1))[end] # measure index: 1-3, bell index: 1-num_bell
res2 = apply_as_qc!(state2, BellMeasure(1,1))[end]

@assert res == res2 # measurement results are the same
