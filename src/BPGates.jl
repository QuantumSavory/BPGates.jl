module BPGates

using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits

export BellState,
    BellSinglePermutation, BellDoublePermutation, BellPauliPermutation,
    BellMeasure, bellmeasure!,
    BellGateQC,
    rand_state,
    apply_as_qc!, convert2QC

function int_to_bit(int,digits)
    int = int - 1 # -1 so that we use julia indexing convenctions
    Bool[int>>shift&0x1 for shift in 0:digits-1]
end
@inline function int_to_bit(int,::Val{2}) # faster if we know we need only two bits
    int = int - 1 # -1 so that we use julia indexing convenctions
    int & 0x1, int>>1 & 0x1
end
@inline function int_to_bit(int,::Val{4}) # faster if we know we need only two bits
    int = int - 1 # -1 so that we use julia indexing convenctions
    int & 0x1, int>>1 & 0x1, int>>2 & 0x1, int>>3 & 0x1
end
function bit_to_int(bits)
    reduce(⊻,(bit<<(index-1) for (index,bit) in enumerate(bits))) + 1 # +1 so that we use julia indexing convenctions
end
# faster version when we know how many bits we have (because we know do not need to iterate)
@inline bit_to_int(bit1,bit2) = bit1 ⊻ bit2<<1 + 1 # +1 so that we use julia indexing convenctions
@inline bit_to_int(bit1,bit2,bit3,bit4) = bit1 ⊻ bit2<<1 ⊻ bit3<<2 ⊻ bit4<<3 + 1 # +1 so that we use julia indexing convenctions

"""
A diagonal representation of bell diagonal states that only tracks the phases in front of the stabilizers tableau
instead of the whole stabilizer tableau.

For example, XX -ZZ is represented as 01.

```jldoctest
julia> new_state = BellState([0,1,1,0])
BellState(Bool[0, 1, 1, 0])

julia> new_state.phases
4-element BitVector:
 0
 1
 1
 0
```
"""
struct BellState
    phases::BitVector
end

Base.copy(state::BellState) = BellState(copy(state.phases))

abstract type BellOp <: QuantumClifford.AbstractCliffordOperator end

"""
One type of bell preserving gates that only performs a Pauli permutation on one side (in this case, we apply the Pauli 
gates to Alice side of Bell pairs) of the shared Bell pairs
pidx can take value of 1-4 representating one of the four pauli permutations available
sidx is a tuple that indicates which pair of Bell states the Pauli permutation will be applied to

```jldoctest
julia> BellPauliPermutation(2,(2,3))
BellPauliPermutation(2, (2, 3))
```
"""
struct BellPauliPermutation <: BellOp
    pidx::Int
    sidx::Tuple{Int,Int}
end

"""
One type of bell preserving gates that only performs a single permutation
pidx can take value of 1-6 representating one of the six single qubit permutation (which is a subset of Clifford gates on 1 qubit)
sidx is an integer of which Bell state the permutation will be applied to. Both Alice and Bob will apply the same one qubit Clifford gates
"""
struct BellSinglePermutation <: BellOp
    pidx::Int
    sidx::Int
end

"""
One type of bell preserving gates that only performs a double permutation
pidx can take value of 1-20 representating one of the six double qubit permutation (which is a subset of Clifford gates on 2 qubit)
sidx is an integer of which pair of Bell states the permutation will be applied to. Both Alice and Bob will apply the same two qubit Clifford gates.
"""
struct BellDoublePermutation <: BellOp
    pidx::Int
    sidx::Tuple{Int,Int}
end

"""
Coincidence measurement to detect errors
midx can take value of 1-3 representing measurement performed in x, y, z basis respectively
sidx is an integer of which Bell state the measurement is applied to. The state will be reset to 00 after being applied measurement.
"""
struct BellMeasure <: QuantumClifford.AbstractMeasurement
    midx::Int
    sidx::Int
end

##############################
# Permutations
##############################

const one_perm_tuple = (
    (1, 2, 3, 4),
    (1, 3, 2, 4),
    (3, 1, 2, 4),
    (3, 2, 1, 4),
    (2, 3, 1, 4),
    (2, 1, 3, 4)
)

function QuantumClifford.apply!(state::BellState, op::BellSinglePermutation)
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx*2-1],phase[op.sidx*2])
    @inbounds perm = one_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[op.sidx*2-1] = bit1
    @inbounds phase[op.sidx*2] = bit2
    return state
end

const pauli_perm_tuple = (
    (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
    (3, 4, 1, 2, 7, 8, 5, 6, 11, 12, 9, 10, 15, 16, 13, 14),
    (2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15),
    (4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 16, 15, 14, 13)
)

function QuantumClifford.apply!(state::BellState, op::BellPauliPermutation)
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx[1]*2-1],phase[op.sidx[1]*2],phase[op.sidx[2]*2-1], phase[op.sidx[2]*2])
    @inbounds perm = pauli_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    changed_phases = int_to_bit(permuted_idx, Val(4))
    @inbounds phase[op.sidx[1]*2-1] = changed_phases[1]
    @inbounds phase[op.sidx[1]*2] = changed_phases[2]
    @inbounds phase[op.sidx[2]*2-1] = changed_phases[3]
    @inbounds phase[op.sidx[2]*2] = changed_phases[4]
    return state
end

const double_perm_tuple = (
    (1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16),
    (1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16),
    (1, 2, 11, 12, 6, 5, 16, 15, 9, 10, 3, 4, 14, 13, 8, 7),
    (1, 3, 10, 12, 7, 5, 16, 14, 9, 11, 2, 4, 15, 13, 8, 6),
    (1, 9, 7, 15, 10, 2, 16, 8, 3, 11, 5, 13, 12, 4, 14, 6),
    (1, 3, 6, 8, 5, 7, 2, 4, 11, 9, 16, 14, 15, 13, 12, 10),
    (9, 11, 6, 8, 5, 7, 10, 12, 3, 1, 16, 14, 15, 13, 4, 2),
    (1, 10, 11, 4, 16, 7, 6, 13, 9, 2, 3, 12, 8, 15, 14, 5),
    (1, 11, 5, 15, 6, 16, 2, 12, 3, 9, 7, 13, 8, 14, 4, 10),
    (1, 11, 7, 13, 16, 6, 10, 4, 3, 9, 5, 15, 14, 8, 12, 2),
    (1, 9, 6, 14, 2, 10, 5, 13, 11, 3, 16, 8, 12, 4, 15, 7),
    (3, 11, 6, 14, 2, 10, 7, 15, 9, 1, 16, 8, 12, 4, 13, 5),
    (1, 7, 2, 8, 5, 3, 6, 4, 10, 16, 9, 15, 14, 12, 13, 11),
    (3, 5, 2, 8, 7, 1, 6, 4, 10, 16, 11, 13, 14, 12, 15, 9),
    (1, 5, 10, 14, 7, 3, 16, 12, 2, 6, 9, 13, 8, 4, 15, 11),
    (3, 7, 10, 14, 5, 1, 16, 12, 2, 6, 11, 15, 8, 4, 13, 9),
    (9, 2, 5, 14, 10, 1, 6, 13, 7, 16, 11, 4, 8, 15, 12, 3),
    (9, 10, 7, 8, 2, 1, 16, 15, 5, 6, 11, 12, 14, 13, 4, 3),
    (9, 7, 6, 12, 5, 11, 10, 8, 16, 2, 3, 13, 4, 14, 15, 1),
    (10, 6, 3, 15, 5, 9, 16, 4, 11, 7, 2, 14, 8, 12, 13, 1)
)

function QuantumClifford.apply!(state::BellState, op::BellDoublePermutation)
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx[1]*2-1],phase[op.sidx[1]*2],phase[op.sidx[2]*2-1], phase[op.sidx[2]*2])
    @inbounds perm = double_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    changed_phases = int_to_bit(permuted_idx, Val(4))
    @inbounds phase[op.sidx[1]*2-1] = changed_phases[1]
    @inbounds phase[op.sidx[1]*2] = changed_phases[2]
    @inbounds phase[op.sidx[2]*2-1] = changed_phases[3]
    @inbounds phase[op.sidx[2]*2] = changed_phases[4]
    return state
end

##############################
# Measurements
##############################

const measure_tuple = (
    (true, true, true),
    (false, false, true),
    (true, false, false),
    (false, true, false)
)

"""
Apply coincidence measurement on a bell state.
Return state, false if an error is detected
The measured state will be reset to 00.

```jldoctest
julia> bellmeasure!(BellState([0,1,1,1]), BellMeasure(2,1))
(BellState(Bool[0, 0, 1, 1]), false)
```

"""
function bellmeasure!(state::BellState, op::BellMeasure)
    phase = state.phases
    result = measure_tuple[bit_to_int(phase[op.sidx*2-1],phase[op.sidx*2])][op.midx]
    phase[op.sidx*2-1:op.sidx*2] .= 0 # reset the measured pair to 00
    return state, result
end

# TODO remove Experimental.NoisyCircuits namespacing when possible
function QuantumClifford.Experimental.NoisyCircuits.applywstatus!(state::BellState, op::BellMeasure)
    state, result = bellmeasure!(state, op)
    state, result ? QuantumClifford.Experimental.NoisyCircuits.continue_stat : QuantumClifford.Experimental.NoisyCircuits.failure_stat
end

##############################
# Full BP gate
##############################

"""
Most general representation of a Bell preserving gate on two qubits.
The general gate consists of a pauli permutation, a double qubit permutation, two single qubit permutations, 
and two indices indicating which pair of Bell states the general bell preserving gate will be applied to.

```jldoctest
julia> BellGateQC(3,15,2,1,4,5)
BellGateQC(3, 15, 2, 1, 4, 5)
```
"""
struct BellGateQC <: BellOp
    pauli::Int
    double::Int
    single1::Int
    single2::Int
    idx1::Int
    idx2::Int
end

function QuantumClifford.apply!(state, g::BellGateQC)
    apply!(state, BellPauliPermutation(g.pauli,(g.idx1,g.idx2)))
    apply!(state, BellDoublePermutation(g.double,(g.idx1,g.idx2)))
    apply!(state, BellSinglePermutation(g.single1,g.idx1))
    apply!(state, BellSinglePermutation(g.single2,g.idx2))
    return state
end

function QuantumClifford.apply!(state, circuit)
    for gate in circuit
        if apply!(state, gate)[end] == false
            return state, false
        end
    end
    return state, true
end

##############################
# Convertion from BP to QC
##############################

stab2qidx(stab)=isone.(stab.phases.÷0x2)

const one_perm_qc = (
    C"X Z",
    C"Z X",
    C"Y X",
    C"X Y",
    C"Z Y",
    C"Y Z"
)

const two_perm_qc = (
    C"X_ _Z Z_ _X",
    C"_X X_ _Z Z_",
    C"XX _X Z_ ZZ",
    C"ZX _X X_ XZ",
    C"XZ X_ _X ZX",
    C"ZZ XX X_ _Z",
    C"ZY XX X_ _Y",
    C"XX _X ZX YY",
    C"_Z X_ XX ZZ",
    C"XZ X_ XX YY",
    C"ZZ XX _X Z_",
    C"YZ XX _X Y_",
    C"Z_ ZX XZ _Z",
    C"Y_ YX XZ _Z",
    C"ZX Z_ _Z XZ",
    C"YX Y_ _Z XZ",
    C"_Y XY ZX Z_",
    C"XY _Y Z_ ZX",
    C"ZY YZ XY _Y",
    C"YX Y_ _Y ZY"
)

const pauli_perm_qc = (P"II",P"XI",P"ZI",P"YI")

"""
Initialize a random Bell diagonal state
Input is the number of shared Bell pairs in the entanglement network.
"""
function rand_state(num_bell)  # TODO this would fail above 32 Bell pairs
    return BellState(int_to_bit(rand(1:4^num_bell),num_bell*2))
end

"""
Convert a Bell perserving gate or state in BPGates.jl representation to the corresponding 
Clifford gate or stabilizer formalism in QuantumClifford.jl representation

```jldoctest
julia> convert2QC(BellSinglePermutation(2,3))
2-element Vector{Tuple{CliffordOperator{Vector{UInt8}, Matrix{UInt64}}, Vector{Int64}}}:
 (X ⟼ + Z
Z ⟼ + X
, [5])
 (X ⟼ + Z
Z ⟼ + X
, [6])

julia> convert2QC(BellState([1,0,0,1]))
- XX__
+ ZZ__
+ __XX
- __ZZ
```
"""
function convert2QC end

function convert2QC(gate::BellSinglePermutation)
    return [(one_perm_qc[gate.pidx], [gate.sidx*2-1]),
    (one_perm_qc[gate.pidx], [gate.sidx*2])]
end
function convert2QC(gate::BellDoublePermutation)
    return [(two_perm_qc[gate.pidx], [gate.sidx[1]*2-1, gate.sidx[2]*2-1]),
    (two_perm_qc[gate.pidx], [gate.sidx[1]*2, gate.sidx[2]*2])]
end
function convert2QC(gate::BellPauliPermutation)
    return [(pauli_perm_qc[gate.pidx], [gate.sidx[1]*2-1, gate.sidx[2]*2-1])]
end

function convert2QC(gate::BellGateQC)
    return [
        (pauli_perm_qc[gate.pauli], [gate.idx1*2-1, gate.idx2*2-1]),
        (two_perm_qc[gate.double], [gate.idx1*2-1, gate.idx2*2-1]),
        (two_perm_qc[gate.double], [gate.idx1*2, gate.idx2*2]),
        (one_perm_qc[gate.single1], [gate.idx1*2-1]),
        (one_perm_qc[gate.single1], [gate.idx1*2]),
        (one_perm_qc[gate.single2], [gate.idx2*2-1]),
        (one_perm_qc[gate.single2], [gate.idx2*2]),
    ]
end

function convert2QC(state::BellState)
    res = bell(length(state.phases)÷2)
    res.phases .= state.phases*0x2
    return res
end

"""
Take in BPGates.jl state and gate and apply the bell preserving gate on a bell diagonal state in terms of QuantumClifford.jl
formalisms and return results in BPGates.jl formalism.

julia> apply_as_qc!(BellState([0,1,1,0]), BellSinglePermutation(2,1))
(BellState(Bool[1, 0, 1, 0]), true)

julia> apply_as_qc!(BellState([0,1,1,0]), BellMeasure(2,1))
(BellState(Bool[0, 0, 1, 0]), false)
"""
function apply_as_qc! end

function apply_as_qc!(state::BellState, gate::BellOp)
    s = convert2QC(state)
    for (g, idx) in convert2QC(gate)
        apply!(s, g, idx)
    end
    new_phases = [project!(s,proj)[end]÷2 for proj in bell(length(state.phases)÷2)]
    state.phases .= new_phases
    return state, true
end

function apply_as_qc!(state::BellState, gate::BellMeasure)
    phases = state.phases[:]
    s = MixedDestabilizer(convert2QC(state))
    if gate.midx == 1
        res = (projectXrand!(s,gate.sidx*2-1)[2]==projectXrand!(s,gate.sidx*2)[2]) ? true : false
    elseif gate.midx == 2
        res = (projectYrand!(s,gate.sidx*2-1)[2]!=projectYrand!(s,gate.sidx*2)[2]) ? true : false
    else
        res = (projectZrand!(s,gate.sidx*2-1)[2]==projectZrand!(s,gate.sidx*2)[2]) ? true : false
    end
    phases[gate.sidx*2-1:gate.sidx*2] .=0
    state.phases .= phases
    return state, res
end

function apply_as_qc!(state::BellState, circuit)
    for op in circuit
        if apply_as_qc!(state,op)[end]==false
            return state, false
        end
    end
    return state, true
end

end # module
