module BPGates

using Permutations

export Int2BitArray, BitArray2Int, stab2qidx, BellState, BellSinglePermutation, BellDoublePermutation, BellPauliPermutation, BellMeasure, BellGate, apply_op!, initialize

Int2BitArray(x,padding)=reverse(BitArray(digits(x, base=2,pad=padding)))
BitArray2Int(x)=UInt8(reverse(x).chunks[1])
stab2qidx(stab)=BitArray(stab.phases.÷0x2)

pauli_perm_list = [
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
[3, 4, 1, 2, 7, 8, 5, 6, 11, 12, 9, 10, 15, 16, 13, 14],
[2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15],
[4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 16, 15, 14, 13]
]
const pauli_perm = [Permutation(p) for p in pauli_perm_list]

test() = 2

one_perm_list = [[1, 2, 3, 4],
 [1, 3, 2, 4],
 [3, 1, 2, 4],
 [3, 2, 1, 4],
 [2, 3, 1, 4],
 [2, 1, 3, 4]]
const one_perm = [Permutation(p) for p in one_perm_list]

two_perm_list =  [[1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16],
 [1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16],
 [1, 2, 11, 12, 6, 5, 16, 15, 9, 10, 3, 4, 14, 13, 8, 7],
 [1, 3, 10, 12, 7, 5, 16, 14, 9, 11, 2, 4, 15, 13, 8, 6],
 [1, 9, 7, 15, 10, 2, 16, 8, 3, 11, 5, 13, 12, 4, 14, 6],
 [1, 3, 6, 8, 5, 7, 2, 4, 11, 9, 16, 14, 15, 13, 12, 10],
 [9, 11, 6, 8, 5, 7, 10, 12, 3, 1, 16, 14, 15, 13, 4, 2],
 [1, 10, 11, 4, 16, 7, 6, 13, 9, 2, 3, 12, 8, 15, 14, 5],
 [1, 11, 5, 15, 6, 16, 2, 12, 3, 9, 7, 13, 8, 14, 4, 10],
 [1, 11, 7, 13, 16, 6, 10, 4, 3, 9, 5, 15, 14, 8, 12, 2],
 [1, 9, 6, 14, 2, 10, 5, 13, 11, 3, 16, 8, 12, 4, 15, 7],
 [3, 11, 6, 14, 2, 10, 7, 15, 9, 1, 16, 8, 12, 4, 13, 5],
 [1, 7, 2, 8, 5, 3, 6, 4, 10, 16, 9, 15, 14, 12, 13, 11],
 [3, 5, 2, 8, 7, 1, 6, 4, 10, 16, 11, 13, 14, 12, 15, 9],
 [1, 5, 10, 14, 7, 3, 16, 12, 2, 6, 9, 13, 8, 4, 15, 11],
 [3, 7, 10, 14, 5, 1, 16, 12, 2, 6, 11, 15, 8, 4, 13, 9],
 [9, 2, 5, 14, 10, 1, 6, 13, 7, 16, 11, 4, 8, 15, 12, 3],
 [9, 10, 7, 8, 2, 1, 16, 15, 5, 6, 11, 12, 14, 13, 4, 3],
 [9, 7, 6, 12, 5, 11, 10, 8, 16, 2, 3, 13, 4, 14, 15, 1],
 [10, 6, 3, 15, 5, 9, 16, 4, 11, 7, 2, 14, 8, 12, 13, 1]]
const two_perm = [Permutation(p) for p in two_perm_list]

abstract type BellOp end

struct BellState
    phases::BitVector
end

struct BellPauliPermutation <:BellOp
    pidx::UInt
    sidx
end

struct BellSinglePermutation <:BellOp
    pidx::UInt
    sidx
end

struct BellDoublePermutation <:BellOp
    pidx::UInt
    sidx
end

struct BellMeasure <: BellOp
    midx::Int
    sidx
end

Base.getindex(a::AbstractVector, b::BitVector)=a[b.chunks[1]]

##############################
# Permutations
##############################

function apply_op!(state::BellState, op::BellSinglePermutation)
    phase = state.phases
    phase[op.sidx*2-1:op.sidx*2] = Int2BitArray(one_perm[op.pidx][BitArray2Int(phase[op.sidx*2-1:op.sidx*2])+1]-1,2)
    return BellState(phase),:continue
end

function apply_op!(state::BellState, op::BellDoublePermutation) 
    phase = state.phases
    changed_phases = Int2BitArray(two_perm[op.pidx][BitArray2Int(vcat(phase[op.sidx[1]*2-1:op.sidx[1]*2],phase[op.sidx[2]*2-1:op.sidx[2]*2]))+1]-1,4)
    phase[op.sidx[1]*2-1:op.sidx[1]*2] = changed_phases[1:2]
    phase[op.sidx[2]*2-1:op.sidx[2]*2] = changed_phases[3:4]
    return BellState(phase),:continue
end

function apply_op!(state::BellState, op::BellPauliPermutation)
    phase = state.phases
    changed_phases = Int2BitArray(pauli_perm[op.pidx][BitArray2Int(vcat(phase[op.sidx[1]*2-1:op.sidx[1]*2],phase[op.sidx[2]*2-1:op.sidx[2]*2]))+1]-1,4)
    phase[op.sidx[1]*2-1:op.sidx[1]*2] = changed_phases[1:2]
    phase[op.sidx[2]*2-1:op.sidx[2]*2] = changed_phases[3:4]
    return BellState(phase),:continue
end

##############################
# Bell Gates
##############################

struct BellGate
    pauli::BellPauliPermutation
    double::BellDoublePermutation
    single1::BellSinglePermutation
    single2::BellSinglePermutation
end

BellGate(pauli, double, single1, single2, idx1, idx2) = BellGate(BellPauliPermutation(pauli,[idx1,idx2]),
    BellDoublePermutation(double,[idx1,idx2]),
    BellSinglePermutation(single1,idx1),
    BellSinglePermutation(single2,idx2))

function apply_op!(state,bell_gate::BellGate)
    circuit = [bell_gate.pauli, bell_gate.double, bell_gate.single1, bell_gate.single2]
    for op in circuit
        if apply_op!(state,op)[end]==:detected_failure
            return state, :detected_failure
        end
    end
    return state, :continue
end

##############################
# Measurements
##############################

mresult=falses(4,3)
for i in 1:3 
    mresult[1,i]=true
    mresult[i+1,i]=true
end

const rmeasurement=mresult

function apply_op!(state::BellState, op::BellMeasure) 
    phase = state.phases
    result=apply_op!(phase[op.sidx*2-1:op.sidx*2], op)
    phase[op.sidx*2-1:op.sidx*2]=BitVector([0,0])
    if result
        return BellState(phase), :continue
    else
        return BellState(phase), :detected_failure
    end
end

function apply_op!(phase, op::BellMeasure)
    result=rmeasurement[BitArray2Int(phase)+1,op.midx]
    return result
end

function apply_op!(state,circuit)
    for op in circuit
        if apply_op!(state,op)[end]==:detected_failure
            return state, :detected_failure
        end
    end
    return state, :continue
end

function initialize(num_bell) 
    return BellState(Int2BitArray(rand(1:4^num_bell),num_bell*2)) 
end

end # module
