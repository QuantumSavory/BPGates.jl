import QuantumClifford: affectedqubits, noisify, skip_idling_noise
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

# plain noise structs to just capture noise info, not ops
struct PauliNoiseData
    px::Float64
    py::Float64
    pz::Float64
end

struct T1NoiseData
    λ₁::Float64
end

struct T2NoiseData
    λ₂::Float64
end
# functions to build noise ops out of noise data, helper function used in idle noise 
build_noise_op(idx::Int, n::PauliNoiseData) = PauliNoiseOp(idx, n.px, n.py, n.pz)
build_noise_op(idx::Int, n::T1NoiseData)    = T1NoiseOp(idx, n.λ₁)
build_noise_op(idx::Int, n::T2NoiseData)    = T2NoiseOp(idx, n.λ₂)
const BPNoiseData = Union{PauliNoiseData, T1NoiseData, T2NoiseData}

# a circuit noise equivalent for bpgates, similar to what is implemented in noisify 

struct BPCircuitNoise
    single_pair::Union{BPNoiseData, Nothing}
    two_pair::Union{BPNoiseData, Nothing}
    idle_noise::Union{Vector{BPNoiseData}, Nothing}
    measurement::Union{Float64, Nothing}   # flip probability p, used by NoisyBellMeasure
end

function BPCircuitNoise(
    single_pair::Union{BPNoiseData, Nothing} = nothing,
    two_pair::Union{BPNoiseData, Nothing}    = nothing,
    idle_noise::Union{Vector{BPNoiseData}, Nothing}  = nothing;
    measurement::Union{Float64, Nothing}      = nothing,
)
    return BPCircuitNoise(single_pair, two_pair, idle_noise, measurement)
end

BPCircuitNoise(noise::BPNoiseData;measurement=nothing) = BPCircuitNoise(
    single_pair = noise,
    two_pair    = noise,
    idle_noise   = [noise],
    measurement  = measurement
)

# for idle noise 
skip_idling_noise(op) = false  
skip_idling_noise(::PauliNoiseOp) = true
skip_idling_noise(::T1NoiseOp) = true
skip_idling_noise(::T2NoiseOp) = true
skip_idling_noise(::PauliNoiseBellGate) = true
skip_idling_noise(::NoisyBellMeasure) = true
skip_idling_noise(::NoisyBellMeasureNoisyReset) = true



insert_idle_noise_bp(circuit::AbstractVector, ::Nothing) = circuit

function insert_idle_noise_bp(circuit::AbstractVector, idle_noise::BPNoiseData)
    isempty(circuit) && return copy(circuit)

    all_qubits = [q for op in circuit
                    if !skip_idling_noise(op)
                    for q in affectedqubits(op)]
    isempty(all_qubits) && return copy(circuit)

    first_seen    = Dict{Int, Int}()
    filled_up_to  = Dict{Int, Int}()
    op_steps      = Int[]
    active_qubits = Dict{Int, Set{Int}}()

    for op in circuit
        if skip_idling_noise(op)
            push!(op_steps, 0)
            continue
        end

        qs = collect(affectedqubits(op))
        if isempty(qs)
            push!(op_steps, 0)
            continue
        end

        for q in qs
            get!(filled_up_to, q, 1)
        end

        step = maximum(filled_up_to[q] for q in qs)
        push!(op_steps, step)
        union!(get!(active_qubits, step, Set{Int}()), qs)

        for q in qs
            get!(first_seen, q, step)
            filled_up_to[q] = step + 1
        end
    end

    output  = []
    emitted = Set{Int}()

    for (op, step) in zip(circuit, op_steps)
        if !skip_idling_noise(op) && !(step in emitted)
            active = get(active_qubits, step, Set{Int}())

            idle = sort([q for q in keys(first_seen) if !(q in active) && first_seen[q] < step])

            for q in idle
                for n in idle_noise
                    push!(output, build_noise_op(q, n))
                end
            end

            push!(emitted, step)
        end

        push!(output, op)
    end

    return output
end



function noisify(circuit::AbstractVector, noise_model::BPCircuitNoise)
    idle_noisy_circuit = insert_idle_noise_bp(circuit, noise_model.idle_noise)
    return reduce(vcat, noisify.(idle_noisy_circuit, (noise_model,)))
end


noisify(op, ::BPCircuitNoise) = [op]
noisify(op, ::BPNoiseData) = [op]
noisify(op, ::Nothing) = [op]

# Already-noisy ops pass through (Very messy since there's no specific type for noisy ops)
noisify(op::PauliNoiseOp,               ::BPCircuitNoise) = [op]
noisify(op::T1NoiseOp,                  ::BPCircuitNoise) = [op]
noisify(op::T2NoiseOp,                  ::BPCircuitNoise) = [op]
noisify(op::PauliNoiseBellGate,         ::BPCircuitNoise) = [op]
noisify(op::NoisyBellMeasure,           ::BPCircuitNoise) = [op]
noisify(op::NoisyBellMeasureNoisyReset, ::BPCircuitNoise) = [op]



# Initial dispatch methods (first methods that noisify goes through)
noisify(op::BellSinglePermutation,  m::BPCircuitNoise) = noisify(op, m.single_pair)
noisify(op::BellPauliPermutation,   m::BPCircuitNoise) = noisify(op, m.single_pair)
noisify(op::BellDoublePermutation,  m::BPCircuitNoise) = noisify(op, m.two_pair)
noisify(op::BellGate,               m::BPCircuitNoise) = noisify(op, m.two_pair)
noisify(op::BellSwap,               m::BPCircuitNoise) = noisify(op, m.two_pair)
noisify(op::BellMeasure,            m::BPCircuitNoise) = noisify(op, m.measurement)
noisify(op::CNOTPerm, m::BPCircuitNoise) = noisify(op, m.two_pair)


# sub dispatch methods, the next set of dispatch methods after the initial methods 
noisify(op::BellSinglePermutation, n::BPNoiseData) = [op, build_noise_op(only(affectedqubits(op)), n)]

noisify(op::BellPauliPermutation, n::BPNoiseData) = [op, build_noise_op(only(affectedqubits(op)), n)]

noisify(op::BellDoublePermutation, n::BPNoiseData) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]

noisify(op::BellGate, n::BPNoiseData) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]

noisify(op::BellSwap, n::BPNoiseData) = [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::CNOTPerm, n::BPNoiseData) =  [op, build_noise_op(affectedqubits(op)[1], n), build_noise_op(affectedqubits(op)[2], n)]
noisify(op::BellMeasure, p::Float64) = [NoisyBellMeasureNoisyReset(op, p)]
noisify(op::BellMeasure, ::Nothing)  = [op]

