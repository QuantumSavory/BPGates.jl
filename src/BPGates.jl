module BPGates

using Permutations
using QuantumClifford

export int_to_bit, bit_to_int, copy_state, BellState, BellSinglePermutation, BellDoublePermutation, BellPauliPermutation, BellMeasure, BellGateQC, apply_op!, rand_state, stab2qidx, apply_as_qc!, convert2QC

function int_to_bit(int,digits)
    int = int - 1 # -1 so that we use julia indexing convenctions
    Bool[int>>shift&0x1 for shift in 0:digits-1]    
end
function bit_to_int(bits)
    reduce(⊻,(bit<<(index-1) for (index,bit) in enumerate(bits))) + 1 # +1 so that we use julia indexing convenctions
end

@assert all(bit_to_int(int_to_bit(n,4))==n for n in 1:2^4)

pauli_perm_list = [
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
[3, 4, 1, 2, 7, 8, 5, 6, 11, 12, 9, 10, 15, 16, 13, 14],
[2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15],
[4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 16, 15, 14, 13]
]
const pauli_perm = [Permutation(p) for p in pauli_perm_list]

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
# Helper functions
##############################

function copy_state(state::BellState)
    return BellState(copy(state.phases))    
end

##############################
# Permutations
##############################

function apply_op!(state::BellState, op::BellSinglePermutation)
    phase = state.phases
    phase[op.sidx*2-1:op.sidx*2] = int_to_bit(one_perm[op.pidx][bit_to_int(phase[op.sidx*2-1:op.sidx*2])],2)
    return BellState(phase),:continue
end

function apply_op!(state::BellState, op::BellDoublePermutation) 
    phase = state.phases
    changed_phases = int_to_bit(two_perm[op.pidx][bit_to_int(vcat(phase[op.sidx[1]*2-1:op.sidx[1]*2],phase[op.sidx[2]*2-1:op.sidx[2]*2]))],4)
    phase[op.sidx[1]*2-1:op.sidx[1]*2] = changed_phases[1:2]
    phase[op.sidx[2]*2-1:op.sidx[2]*2] = changed_phases[3:4]
    return BellState(phase),:continue
end

function apply_op!(state::BellState, op::BellPauliPermutation)
    phase = state.phases
    changed_phases = int_to_bit(pauli_perm[op.pidx][bit_to_int(vcat(phase[op.sidx[1]*2-1:op.sidx[1]*2],phase[op.sidx[2]*2-1:op.sidx[2]*2]))],4)
    phase[op.sidx[1]*2-1:op.sidx[1]*2] = changed_phases[1:2]
    phase[op.sidx[2]*2-1:op.sidx[2]*2] = changed_phases[3:4]
    return BellState(phase),:continue
end

##############################
# Measurements
##############################

mresult = falses(4,3)
for idx in [[1,1],[1,2],[1,3],[2,3],[3,1],[4,2]]
    mresult[idx[1],idx[2]] = true
end

const rmeasurement=mresult

function apply_op!(state::BellState, op::BellMeasure) 
    phase = state.phases
    result = apply_op!(phase[op.sidx*2-1:op.sidx*2], op)
    phase[op.sidx*2-1:op.sidx*2]=BitVector([0,0]) # reset the measured pair to 00
    if result
        return BellState(phase), :continue
    else
        return BellState(phase), :detected_error
    end
end

function apply_op!(phase, op::BellMeasure)
    result=rmeasurement[bit_to_int(phase),op.midx]
    return result
end

##############################
# Full BP gate
##############################

struct BellGateQC
    pauli
    double
    single1
    single2
    idx1
    idx2
end

function apply_op!(state, g::BellGateQC)
    circuit = [
    BellPauliPermutation(g.pauli,[g.idx1,g.idx2]),
    BellDoublePermutation(g.double,[g.idx1,g.idx2]),
    BellSinglePermutation(g.single1,g.idx1),
    BellSinglePermutation(g.single2,g.idx2)
    ]
    for op in circuit
        if apply_op!(state,op)[end]==:detected_error
            return state, :detected_error
        end
    end
    return state, :continue
end

function apply_op!(state,circuit)
    for op in circuit
        if apply_op!(state,op)[end]==:detected_error
            return state, :detected_error
        end
    end
    return state, :continue
end

##############################
# Convertion from BP to QC
##############################

stab2qidx(stab)=isone.(stab.phases.÷0x2)

const one_perm_qc = [C"X Z",
C"Z X",
C"Y X",
C"X Y",
C"Z Y",
C"Y Z"];

const two_perm_qc = [C"X_ _Z Z_ _X",
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
C"YX Y_ _Y ZY"];

const pauli_perm_qc = [P"II",P"XI",P"ZI",P"YI"];

function rand_state(num_bell) 
    return BellState(int_to_bit(rand(1:4^num_bell),num_bell*2)) 
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
function convert2BP(state::Stabilizer)
    return BellState(stab2qidx(state))
end

function apply_as_qc!(state::BellState, gate::BellGateQC)
    s = convert2QC(state)
    for (g, idx) in convert2QC(gate)
        apply!(s, g, idx)
    end
    new_phases = [project!(s,proj)[end]÷2 for proj in bell(length(state.phases)÷2)]
    state.phases .= new_phases
    return state, :continue
end


function apply_as_qc!(state::BellState, gate::BellMeasure)
    phases = state.phases[:]
    s = MixedDestabilizer(convert2QC(state))
    if gate.midx == 1
        res = (projectXrand!(s,gate.sidx*2-1)[2]==projectXrand!(s,gate.sidx*2)[2]) ? :continue : :detected_error
        
    elseif gate.midx == 2
        res = !(projectYrand!(s,gate.sidx*2-1)[2]==projectYrand!(s,gate.sidx*2)[2]) ? :continue : :detected_error
    else
        res = (projectZrand!(s,gate.sidx*2-1)[2]==projectZrand!(s,gate.sidx*2)[2]) ? :continue : :detected_error
    end
    phases[gate.sidx*2-1:gate.sidx*2] .=0
    state.phases .= phases
    return state, res
end

function apply_as_qc!(state::BellState, circuit)
    for op in circuit
        if apply_as_qc!(state,op)[end]==:detected_error
            return state, :detected_error
        end
    end
    return state, :continue
end

end # module
