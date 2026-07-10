import QuantumClifford: affectedqubits, noisify, skip_idling_noise, append_idle_noise!, AbstractNoise, PauliNoise
using QuantumClifford: insert_idle_noise
# to separate measurement ops from other types

# get affectedqubits information to feed into idle noise 
affectedqubits(op::BellSinglePermutation)  = (op.sidx,)
affectedqubits(op::BellPauliPermutation)   = (op.sidx,)
affectedqubits(op::BellDoublePermutation)  = (op.sidx1, op.sidx2)
affectedqubits(op::BellGate)               = (op.idx1, op.idx2)
affectedqubits(op::BellSwap)               = (op.idx1, op.idx2)
affectedqubits(op::BellMeasure)            = (op.sidx,)
affectedqubits(op::PauliNoiseOp)           = (op.idx,)
affectedqubits(op::PauliNoiseBellGate)     = (op.g.idx1, op.g.idx2)
affectedqubits(op::NoisyBellMeasure)       = (op.m.sidx,)
affectedqubits(op::NoisyBellMeasureNoisyReset) = (op.m.sidx,)
affectedqubits(op::T1NoiseOp) = (op.idx,)
affectedqubits(op::T2NoiseOp) = (op.idx,)
affectedqubits(op::CNOTPerm) = (op.idx1, op.idx2)
affectedqubits(gate::GoodSingleQubitPerm) = (gate.idx,)

# plain noise structs to just capture noise info, not ops
struct T1Noise <: AbstractNoise
    λ₁::Float64
end
struct T2Noise <: AbstractNoise
    λ₂::Float64
end
# functions to build noise ops out of noise data, helper function used in idle noise 
build_noise_op(idx::Int, n::PauliNoise) = PauliNoiseOp(idx, n.px, n.py, n.pz)

build_noise_op(idx::Int, n::T1Noise) = T1NoiseOp(idx, n.λ₁)

build_noise_op(idx::Int, n::T2Noise) = T2NoiseOp(idx, n.λ₂)
build_noise_op(::Int, n::AbstractNoise) = throw(ArgumentError("BPGates does not support noise type $(typeof(n))."))

# a circuit noise equivalent for bpgates, similar to what is implemented in noisify for QuantumClifford

struct BPCircuitNoise
    gate_noise::Union{AbstractNoise,Nothing}
    idle_noise::Union{AbstractNoise,Nothing}
    measurement::Union{Float64,Nothing} # flip probability
end

function BPCircuitNoise(;
    gate_noise::Union{AbstractNoise,Nothing} = nothing,
    idle_noise::Union{AbstractNoise,Nothing} = nothing,
    measurement::Union{Float64,Nothing} = nothing,
)
    return BPCircuitNoise(
        gate_noise,
        idle_noise,
        measurement,
    )
end

function BPCircuitNoise(
    noise::AbstractNoise;
    measurement::Union{Float64,Nothing} = nothing,
)
    return BPCircuitNoise(
        gate_noise = noise,
        idle_noise = noise,
        measurement = measurement,
    )
end

skip_idling_noise(::PauliNoiseOp) = true
skip_idling_noise(::T1NoiseOp) = true
skip_idling_noise(::T2NoiseOp) = true
skip_idling_noise(::PauliNoiseBellGate) = true
skip_idling_noise(::NoisyBellMeasure) = true
skip_idling_noise(::NoisyBellMeasureNoisyReset) = true

append_idle_noise!(output, q::Int, idle_noise::Union{T1Noise, T2Noise, PauliNoise}) = push!(output, build_noise_op(q, idle_noise))


noisify(op, ::BPCircuitNoise) = [op]

function noisify(circuit::AbstractVector, noise_model::BPCircuitNoise)
    idle_noisy_circuit = insert_idle_noise(circuit, noise_model.idle_noise)
    return reduce(vcat, noisify.(idle_noisy_circuit, (noise_model,)))
end


# Already-noisy ops pass through (Very messy since there's no specific type for noisy ops)
noisify(op::PauliNoiseOp,               ::BPCircuitNoise) = [op]
noisify(op::T1NoiseOp,                  ::BPCircuitNoise) = [op]
noisify(op::T2NoiseOp,                  ::BPCircuitNoise) = [op]
noisify(op::PauliNoiseBellGate,         ::BPCircuitNoise) = [op]
noisify(op::NoisyBellMeasure,           ::BPCircuitNoise) = [op]
noisify(op::NoisyBellMeasureNoisyReset, ::BPCircuitNoise) = [op]



# Initial dispatch methods (first methods that noisify goes through)
noisify(op::BellSinglePermutation,  m::BPCircuitNoise) = noisify(op, m.gate_noise)
noisify(op::BellPauliPermutation,   m::BPCircuitNoise) = noisify(op, m.gate_noise)
noisify(op::BellDoublePermutation,  m::BPCircuitNoise) = noisify(op, m.gate_noise)
noisify(op::BellGate,               m::BPCircuitNoise) = noisify(op, m.gate_noise)
noisify(op::BellSwap,               m::BPCircuitNoise) = noisify(op, m.gate_noise)
noisify(op::CNOTPerm, m::BPCircuitNoise)               = noisify(op, m.gate_noise)
noisify(op::BellMeasure,            m::BPCircuitNoise) = noisify(op, m.measurement)



# sub dispatch methods, the next set of dispatch methods after the initial methods 
noisify(op::BellSinglePermutation, n::AbstractNoise) = [op, build_noise_op(only(affectedqubits(op)), n)]
noisify(op::BellPauliPermutation, n::AbstractNoise) = [op, build_noise_op(only(affectedqubits(op)), n)]
noisify(op::BellDoublePermutation, n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellGate, n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellSwap, n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::CNOTPerm, n::AbstractNoise) =  [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellMeasure, p::Float64) = [NoisyBellMeasureNoisyReset(op, p, 0.0, 0.0, 0.0)]
noisify(op::BellMeasure, ::Nothing)  = [op]

