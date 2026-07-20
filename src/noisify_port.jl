import QuantumClifford: affectedqubits, noisify, skip_idling_noise, append_idle_noise!, AbstractNoise, PauliNoise, CircuitNoise
using QuantumClifford: insert_idle_noise

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
""" T1Noise(λ₁) Amplitude damping (T1 relaxation) with decay probability `λ₁`. Applies bilaterally to both qubits of a Bell pair — see [`T1NoiseOp`](@ref) """
struct T1Noise <: AbstractNoise
    λ₁::Float64
end

""" T2Noise(λ₂) Dephasing (T2 decoherence) with probability `λ₂`. Applies bilaterally t both qubits of a Bell pair — see [`T2NoiseOp`](@ref) """
struct T2Noise <: AbstractNoise
    λ₂::Float64
end

""" Probability `p` that a Bell measurement outcome is flipped. A scalar wrapper so measurement noise can satisfy `CircuitNoise`'s `AbstractNoise` field. """
struct MeasurementFlipNoise <: AbstractNoise
    p::Float64
end
# functions to build noise ops out of noise data, helper function used in idle noise 
build_noise_op(idx::Int, n::PauliNoise) = PauliNoiseOp(idx, n.px, n.py, n.pz)
build_noise_op(idx::Int, n::UnbiasedUncorrelatedNoise) = PauliNoiseOp(idx, n.p/3, n.p/3, n.p/3)
build_noise_op(idx::Int, n::T1Noise) = T1NoiseOp(idx, n.λ₁)
build_noise_op(idx::Int, n::T2Noise) = T2NoiseOp(idx, n.λ₂)
build_noise_op(::Int, n::AbstractNoise) = throw(ArgumentError("BPGates does not support noise type $(typeof(n))."))


# helpers used by insert_idle_noise in QuantumClifford
skip_idling_noise(::PauliNoiseOp) = true
skip_idling_noise(::T1NoiseOp) = true
skip_idling_noise(::T2NoiseOp) = true
skip_idling_noise(::PauliNoiseBellGate) = true
skip_idling_noise(::NoisyBellMeasure) = true
skip_idling_noise(::NoisyBellMeasureNoisyReset) = true

append_idle_noise!(output, q::Int, idle_noise::Union{T1Noise, T2Noise, PauliNoise}) = push!(output, build_noise_op(q, idle_noise))

noisify(op::BellMeasure, n::MeasurementFlipNoise) = noisify(op, n.p)


# "already noisy" pass-throughs 
noisify(op::PauliNoiseOp,       ::CircuitNoise) = [op]
noisify(op::T1NoiseOp,          ::CircuitNoise) = [op]
noisify(op::T2NoiseOp,          ::CircuitNoise) = [op]
noisify(op::PauliNoiseBellGate, ::CircuitNoise) = [op]

function noisify(op::NoisyBellMeasureNoisyReset, m::CircuitNoise)
    isnothing(m.measurement) && return [op]
    measurement_noise = m.measurement
    measurement_noise isa MeasurementFlipNoise ||
        throw(ArgumentError(
            "BPGates measurement noise must be of type MeasurementFlipNoise or nothing"
        ))
    p = op.p + m.measurement.p - 2 * op.p * m.measurement.p
    return [NoisyBellMeasureNoisyReset(op.m, p, op.px, op.py, op.pz)]
end

function noisify(op::NoisyBellMeasure, noise::CircuitNoise)
    isnothing(noise.measurement) && return [op]

    p2 = noise.measurement.p
    p = op.p + p2 - 2 * op.p * p2

    return [NoisyBellMeasure(op.m, p)]
end

# initial dispatch
noisify(op::BellSinglePermutation, m::CircuitNoise) = noisify(op, m.single_qubit)
noisify(op::BellPauliPermutation,  m::CircuitNoise) = noisify(op, m.single_qubit)
noisify(op::BellDoublePermutation, m::CircuitNoise) = noisify(op, m.two_qubit)
noisify(op::BellGate,              m::CircuitNoise) = noisify(op, m.two_qubit)
noisify(op::BellSwap,              m::CircuitNoise) = noisify(op, m.two_qubit)
noisify(op::CNOTPerm,              m::CircuitNoise) = noisify(op, m.two_qubit)
function noisify(op::BellMeasure, m::CircuitNoise)
    isnothing(m.measurement) && return [op]
    return noisify(op, m.measurement.p)
end



# sub dispatch methods, the next set of dispatch methods after the initial methods 
noisify(op::BellSinglePermutation, n::AbstractNoise) = [op, build_noise_op(only(affectedqubits(op)), n)]
noisify(op::BellPauliPermutation,  n::AbstractNoise) = [op, build_noise_op(only(affectedqubits(op)), n)]
noisify(op::BellDoublePermutation, n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellGate,              n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellSwap,              n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::CNOTPerm,              n::AbstractNoise) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellMeasure, p::Float64) = [NoisyBellMeasureNoisyReset(op, p, 0.0, 0.0, 0.0)]
noisify(op::BellMeasure, ::Nothing)  = [op]