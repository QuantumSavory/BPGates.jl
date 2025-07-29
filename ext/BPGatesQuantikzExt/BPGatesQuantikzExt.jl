module BPGatesQuantikzExt

import Quantikz
using BPGates


function small_vert_matrix(numbers, suffix="")
    return "{\\tiny\\begin{matrix}"*join(numbers, "\\\\")*"\\end{matrix}}"*suffix
end
function small_tuples(list_of_tuples)
    return "{\\tiny\\begin{matrix}"*join(join.(list_of_tuples), "\\\\")*"\\end{matrix}}"
end

Quantikz.QuantikzOp(op::BellSinglePermutation) = Quantikz.U(small_vert_matrix(BPGates.one_perm_tuple[op.pidx], "\\ \\mathcal{C}_1^*"), op.sidx)
Quantikz.QuantikzOp(op::BellDoublePermutation) = Quantikz.MultiControlU(
    small_vert_matrix(BPGates.double_perm_tuple[op.pidx][1:4])*"\\ "*
    small_vert_matrix(BPGates.double_perm_tuple[op.pidx][5:8])*"\\ "*
    small_vert_matrix(BPGates.double_perm_tuple[op.pidx][9:12])*"\\ "*
    small_vert_matrix(BPGates.double_perm_tuple[op.pidx][13:16], "\\ \\mathcal{C}_2^*"),
    [op.sidx1, op.sidx2])
Quantikz.QuantikzOp(op::BellPauliPermutation) = Quantikz.U(["I", "X", "Z", "Y"][op.pidx], op.sidx)
Quantikz.QuantikzOp(op::BellMeasure) = Quantikz.Measurement(small_tuples([BPGates.int_to_bit(phase, Val(2)) for phase in 1:4 if BPGates.measure_tuple[phase][op.midx]]),op.sidx)
Quantikz.QuantikzOp(op::BellGate) = Quantikz.MultiControlU("BP_{$(op.pauli1),$(op.pauli2),$(op.double),$(op.single1),$(op.single2)}",[op.idx1, op.idx2])
Quantikz.QuantikzOp(op::CNOTPerm) = Quantikz.MultiControlU(
    small_vert_matrix(BPGates.good_perm_tuple[op.single1])*"\\ "*
    small_vert_matrix(BPGates.good_perm_tuple[op.single2])*"\\ X",
    op.idx1, [op.idx2])
Quantikz.QuantikzOp(op::GoodSingleQubitPerm) = Quantikz.U(small_vert_matrix(BPGates.good_perm_tuple[op.single]), op.idx)

Quantikz.QuantikzOp(op::PauliNoiseOp) = Quantikz.Noise([op.idx])

Quantikz.QuantikzOp(op::PauliNoiseBellGate) = Quantikz.QuantikzOp(op.g)
Quantikz.QuantikzOp(op::NoisyBellMeasure) = Quantikz.QuantikzOp(op.m)
Quantikz.QuantikzOp(op::NoisyBellMeasureNoisyReset) = Quantikz.QuantikzOp(op.m)

# Add t1 and t2 noise op visualizations, to show as a gate with 'T1' and 'T2' labels, and also the λ values in the label
Quantikz.QuantikzOp(op::T1NoiseOp) = Quantikz.U(small_vert_matrix(("\\mathcal{\\lambda}_1",round(op.λ₁, digits=2)), "\\ T_1"), op.idx)
Quantikz.QuantikzOp(op::T2NoiseOp) = Quantikz.U(small_vert_matrix(("\\mathcal{\\lambda}_2",round(op.λ₂, digits=2)), "\\ T_2"), op.idx)

# extend Quantikz.affectedqubits to include the new gates
# Get the qubit pairs involved in an operation
Quantikz.affectedqubits(gate::PauliNoiseBellGate) = [gate.g.idx1, gate.g.idx2]
Quantikz.affectedqubits(gate::NoisyBellMeasureNoisyReset) = [gate.m.sidx,]
Quantikz.affectedqubits(gate::BellMeasure) = [gate.sidx,]
Quantikz.affectedqubits(gate::CNOTPerm) = [gate.idx1, gate.idx2]
Quantikz.affectedqubits(gate::GoodSingleQubitPerm) = [gate.idx,]
Quantikz.affectedqubits(gate::BellSinglePermutation) = [gate.sidx,]
Quantikz.affectedqubits(gate::BellPauliPermutation) = [gate.sidx,]
Quantikz.affectedqubits(gate::T1NoiseOp) = [gate.idx,]
Quantikz.affectedqubits(gate::T2NoiseOp) = [gate.idx,]



end
