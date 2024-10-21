module BPGates

using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits

using Random

export BellState,
    BellSinglePermutation, BellDoublePermutation, BellPauliPermutation,
    BellMeasure, bellmeasure!,
    BellGate, CNOTPerm, GoodSingleQubitPerm,
    PauliNoiseOp, PauliNoiseBellGate, NoisyBellMeasure, NoisyBellMeasureNoisyReset,
    BellSwap, NoisyBellSwap

const IT = Union{Int8,Int16,Int32,Int64,UInt8,UInt16,UInt32,UInt64}

function int_to_bit(int::IT,digits)
    int = int - 1 # -1 so that we use julia indexing conventions
    Bool[int>>shift&0x1 for shift in 0:digits-1]
end
@inline function int_to_bit(int::IT,::Val{2}) # faster if we know we need only two bits
    int = int - 1 # -1 so that we use julia indexing conventions
    int & 0x1, int>>1 & 0x1
end
@inline function int_to_bit(int::IT,::Val{4}) # faster if we know we need only two bits
    int = int - 1 # -1 so that we use julia indexing conventions
    int & 0x1, int>>1 & 0x1, int>>2 & 0x1, int>>3 & 0x1
end
function bit_to_int(bits)
    reduce(⊻,(bit<<(index-1) for (index,bit) in enumerate(bits))) + 1 # +1 so that we use julia indexing convenctions
end
# faster version when we know how many bits we have (because we know do not need to iterate)
@inline bit_to_int(bit1,bit2) = bit1 ⊻ bit2<<1 + 1 # +1 so that we use julia indexing conventions
@inline bit_to_int(bit1,bit2,bit3,bit4) = bit1 ⊻ bit2<<1 ⊻ bit3<<2 ⊻ bit4<<3 + 1 # +1 so that we use julia indexing convenctions

"""
A diagonal representation of Bell diagonal states that only tracks the phases in front of
the stabilizers tableau instead of the whole stabilizer tableau.

Capable of representing only tensor products of one or more of the states
```
±XX
±ZZ
```
by tracking only the phases in the first column.
For example, `XX -ZZ` is represented as the bitstring `01`.

This representation permits drastically
faster simulation of entanglement purification circuits.

The `BellState(n)` constructor will create `n` Bell pairs.

```jldoctest
julia> bell_state = BellState([0,1,1,0])
BellState(Bool[0, 1, 1, 0])

julia> Stabilizer(bell_state)
+ XX__
- ZZ__
- __XX
+ __ZZ
```

As mentioned above, we can represent only Bell states. Here is the basis being used:

| `BPGates` notation| Stabilizer tableaux | Kets | in X basis | in Y basis |
|:---|:---|:---|:---|:---|
|`00`|`+XX +ZZ`|`∣00⟩+∣11⟩`|`∣++⟩+∣--⟩`|`∣i₊i₋⟩+∣i₋i₊⟩`|
|`01`|`+XX -ZZ`|`∣01⟩+∣10⟩`|`∣++⟩-∣--⟩`|`∣i₊i₊⟩-∣i₋i₋⟩`|
|`10`|`-XX +ZZ`|`∣00⟩-∣11⟩`|`∣+-⟩+∣-+⟩`|`∣i₊i₊⟩+∣i₋i₋⟩`|
|`11`|`-XX -ZZ`|`∣01⟩-∣10⟩`|`∣+-⟩-∣-+⟩`|`∣i₊i₋⟩-∣i₋i₊⟩`|

You can convert between these descriptions using
- `BPGates` to stabilizer state with `QuantumClifford.Stabilizer(bpgates_state)`
- stabilizer state to ket with `QuantumOptics.Ket`
"""
struct BellState <: QuantumClifford.AbstractStabilizer
    phases::BitVector
end
BellState(n::Integer) = BellState(BitVector(falses(2n)))
BellState(t::Tuple) = BellState(BitVector(t))

Base.copy(state::BellState) = BellState(copy(state.phases))
Base.:(==)(l::BellState,r::BellState) = l.phases==r.phases

abstract type BellOp <: QuantumClifford.AbstractCliffordOperator end

"""
Bell preserving gate performing one of 4 possible "Pauli permutations" on a single Bell pair.

Equivalent to applying a Pauli gate to Alice's side of the Bell pair.

The first argument, `pidx`, specifies the permutation (between 1 and 4).
The second argument, `sidx` indicates which Bell pair is acted on.

```jldoctest
julia> BellPauliPermutation(1,1)*BellState(1) |> Stabilizer
+ XX
+ ZZ

julia> BellPauliPermutation(2,1)*BellState(1) |> Stabilizer
+ XX
- ZZ

julia> BellPauliPermutation(3,1)*BellState(1) |> Stabilizer
- XX
+ ZZ

julia> BellPauliPermutation(4,1)*BellState(1) |> Stabilizer
- XX
- ZZ
```
"""
struct BellPauliPermutation <: BellOp
    pidx::Int
    sidx::Int
    function BellPauliPermutation(p,s)
        1 <= p <= 4 || throw(ArgumentError("The permutation index needs to be between 1 and 4"))
        s>0 || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        new(p,s)
    end
end

"""
Bell preserving gate performing one of 6 possible "Clifford phaseless permutations" on a single Bell pair.

Equivalent to applying certain single-qubit Clifford gates to both Alice and Bob.

The first argument, `pidx`, specifies the permutation (between 1 and 6).
The second argument, `sidx` indicates which Bell pair is acted on.
"""
struct BellSinglePermutation <: BellOp
    pidx::Int
    sidx::Int
    function BellSinglePermutation(p,s)
        1 <= p <= 6 || throw(ArgumentError("The permutation index needs to be between 1 and 6"))
        s>0 || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        new(p,s)
    end
end

"""
Bell preserving gate performing one of 20 possible "Clifford phaseless two-pair permutations" on a two Bell pairs.

Equivalent to applying the same two-qubit Clifford gate to both Alice's and Bob's half-pairs.

The first argument, `pidx`, specifies the permutation (between 1 and 20).
The second argument, `sidx` indicates which Bell pairs are acted on.
"""
struct BellDoublePermutation <: BellOp
    pidx::Int
    sidx1::Int
    sidx2::Int
    function BellDoublePermutation(p,s1,s2)
        1 <= p <= 20 || throw(ArgumentError("The permutation index needs to be between 1 and 20"))
        (s1>0 && s2>0) || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        new(p,s1,s2)
    end
end

"""
Coincidence measurement on Bell pairs.

The first argument, `midx`, specifies the X, Y, Z basis respectively.
The second argument, `sidx`, indicates which Bell pair is being measured.

The state will be reset to `00` after being applied measurement.
"""
struct BellMeasure <: QuantumClifford.AbstractMeasurement
    midx::Int
    sidx::Int
    function BellMeasure(p,s)
        1 <= p <= 3 || throw(ArgumentError("The basis measurement index needs to be between 1 and 3"))
        new(p,s)
    end
end

function Base.:(*)(op::BellOp, s::BellState; phases::Bool=true)
    s = copy(s)
    apply!(s,op)
end

##############################
# Permutations
##############################

"""The permutations realized by [`BellSinglePermutation`](@ref)."""
const one_perm_tuple = (
    (1, 2, 3, 4), # good 1
    (1, 3, 2, 4), # good 3
    (3, 1, 2, 4),
    (3, 2, 1, 4),
    (2, 3, 1, 4),
    (2, 1, 3, 4)
)

"""The permutations realized by [`BellSinglePermutation`](@ref) as Clifford operations."""
const one_perm_qc = ( # TODO switch to symbolic gates
    C"X Z",
    C"Z X",
    C"Y X",
    C"X Y",
    C"Z Y",
    C"Y Z"
)

function QuantumClifford.apply!(state::BellState, op::BellSinglePermutation) # TODO abstract away the permutation application as it is used by other gates too
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx*2-1],phase[op.sidx*2])
    @inbounds perm = one_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[op.sidx*2-1] = bit1
    @inbounds phase[op.sidx*2] = bit2
    return state
end

"""The permutations realized by [`BellPauliPermutation`](@ref)."""
const pauli_perm_tuple = (
    (1, 2, 3, 4),
    (3, 4, 1, 2), ## X flip
    (2, 1, 4, 3), ## Z flip
    (4, 3, 2, 1)  ## Y flip
)

"""The permutations realized by [`BellPauliPermutation`](@ref) as Pauli operations."""
const pauli_perm_qc = (sId1,sX,sZ,sY)

function QuantumClifford.apply!(state::BellState, op::BellPauliPermutation) # TODO abstract away the permutation application as it is used by other gates too
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx*2-1],phase[op.sidx*2])
    @inbounds perm = pauli_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[op.sidx*2-1] = bit1
    @inbounds phase[op.sidx*2] = bit2
    return state
end

"""The permutations realized by [`BellDoublePermutation`](@ref)."""
const double_perm_tuple = (
    (1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16),
    (1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16),
    (1, 2, 11, 12, 6, 5, 16, 15, 9, 10, 3, 4, 14, 13, 8, 7), # bilateral CNOT
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

"""The permutations realized by [`BellDoublePermutation`](@ref) as Clifford operations."""
const double_perm_qc = ( # TODO switch to symbolic gates
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

function QuantumClifford.apply!(state::BellState, op::BellDoublePermutation) # TODO abstract away the permutation application as it is used by other gates too
    phase = state.phases
    @inbounds phase_idx = bit_to_int(phase[op.sidx1*2-1],phase[op.sidx1*2],phase[op.sidx2*2-1], phase[op.sidx2*2])
    @inbounds perm = double_perm_tuple[op.pidx]
    @inbounds permuted_idx = perm[phase_idx]
    changed_phases = int_to_bit(permuted_idx, Val(4))
    @inbounds phase[op.sidx1*2-1] = changed_phases[1]
    @inbounds phase[op.sidx1*2] = changed_phases[2]
    @inbounds phase[op.sidx2*2-1] = changed_phases[3]
    @inbounds phase[op.sidx2*2] = changed_phases[4]
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
Apply coincidence measurement on a Bell state.

Return state and measurement result (`false` if an error is detected).

The measured state will be reset to 00.

```jldoctest
julia> bellmeasure!(BellState([0,1,1,1]), BellMeasure(2,1))
(BellState(Bool[0, 0, 1, 1]), false)
```
"""
function bellmeasure!(state::BellState, op::BellMeasure) # TODO document which index corresponds to which measurement
    phase = state.phases
    result = measure_tuple[bit_to_int(phase[op.sidx*2-1],phase[op.sidx*2])][op.midx]
    phase[op.sidx*2-1:op.sidx*2] .= 0 # reset the measured pair to 00
    return state, result
end

function QuantumClifford.applywstatus!(state::BellState, op::BellMeasure)
    state, result = bellmeasure!(state, op)
    state, result ? continue_stat : failure_stat
end

##############################
# Full BP gate
##############################

"""
Most general representation of a Bell preserving gate on two qubits.
The general gate consists of a two Pauli permutations, a double qubit permutation,
two single qubit permutations, and two indices indicating which pair of Bell states the
general Bell preserving gate will be applied to.
"""
struct BellGate <: BellOp
    pauli1::Int
    pauli2::Int
    double::Int
    single1::Int
    single2::Int
    idx1::Int
    idx2::Int
    function BellGate(p1,p2,d,s1,s2,i1,i2)
        (1 <= p1 <= 4 && 1 <= p2 <= 4) || throw(ArgumentError("The Pauli permutation index needs to be between 1 and 4."))
        (1 <= s1 <= 6 && 1 <= s2 <= 6) || throw(ArgumentError("The single-qubit Clifford permutation index needs to be between 1 and 6."))
        (1 <= d <= 20) || throw(ArgumentError("The double-pair Clifford permutation index needs to be between 1 and 20."))
        (i1 > 0 && i2 > 0) || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        i1 != i2 || throw(ArgumentError("The gate has to act on two different Bell pairs, i.e. idx1!=idx2."))
        new(p1,p2,d,s1,s2,i1,i2)
    end
end

function QuantumClifford.apply!(state::BellState, g::BellGate)
    apply!(state, BellPauliPermutation(g.pauli1,g.idx1))
    apply!(state, BellPauliPermutation(g.pauli2,g.idx2))
    apply!(state, BellDoublePermutation(g.double,g.idx1,g.idx2))
    apply!(state, BellSinglePermutation(g.single1,g.idx1))
    apply!(state, BellSinglePermutation(g.single2,g.idx2))
    return state
end

##############################
# SWAP gate
##############################

"""SWAP gate"""
struct BellSwap <: BellOp
    idx1::Int
    idx2::Int
end

function QuantumClifford.apply!(state::BellState, op::BellSwap)
    phase = state.phases
    @inbounds temp1 = phase[op.idx1*2-1]
    @inbounds temp2 = phase[op.idx1*2]
    @inbounds phase[op.idx1*2-1] = phase[op.idx2*2-1]
    @inbounds phase[op.idx1*2] = phase[op.idx2*2]
    @inbounds phase[op.idx2*2-1] = temp1
    @inbounds phase[op.idx2*2] =temp2
    return state
end

##############################
# Typically good operations
##############################

"""A bilateral CNOT preceded by permutations on each of the pairs that map the `00` state to itself."""
struct CNOTPerm <: BellOp
    single1::Int
    single2::Int
    idx1::Int
    idx2::Int
    function CNOTPerm(s1,s2,i1,i2)
        (1 <= s1 <= 6 && 1 <= s2 <= 6) || throw(ArgumentError("The permutation index needs to be between 1 and 6."))
        (i1 > 0 && i2 > 0) || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        i1 != i2 || throw(ArgumentError("The gate has to act on two different Bell pairs, i.e. idx1!=idx2."))
        new(s1,s2,i1,i2)
    end
end

"""A single-qubit Clifford operation acting as a permutation that maps the `00` state to itself."""
struct GoodSingleQubitPerm <: BellOp
    single::Int
    idx::Int
    function GoodSingleQubitPerm(s,i)
        (1 <= s <= 6) || throw(ArgumentError("The permutation index needs to be between 1 and 6."))
        (i > 0) || throw(ArgumentError("The Bell pair index has to be a positive integer."))
        new(s,i)
    end
end

const good_perm_tuple = (
    (1, 2, 3, 4), # perm 1
    (1, 2, 4, 3), #
    (1, 3, 2, 4), # perm 2
    (1, 3, 4, 2), #
    (1, 4, 2, 3), #
    (1, 4, 3, 2), #
)

const h = tHadamard
const p = tPhase
const hp = h*p
const ph = p*h
const good_perm_qc = ( # From the appendix of Optimized Entanglement Purification, but be careful with index notation being different
    (tId1,tId1), # TODO switch to symbolic gates
    (h*ph*ph,h*hp*hp*hp*hp),
    (h,h),
    (ph*ph,hp*hp*hp*hp),
    (ph,hp*hp),
    (h*ph,p*hp)
)

const cnot_perm = (1, 2, 11, 12, 6, 5, 16, 15, 9, 10, 3, 4, 14, 13, 8, 7)

function QuantumClifford.apply!(state::BellState, g::CNOTPerm) # TODO abstract away the permutation application as it is used by other gates too
    phase = state.phases
    @inbounds phase_idxa = bit_to_int(phase[g.idx1*2-1],phase[g.idx1*2])
    @inbounds phase_idxb = bit_to_int(phase[g.idx2*2-1],phase[g.idx2*2])
    if phase_idxa==phase_idxb==1
        return state
    end
    # first qubit permutation
    @inbounds perm = good_perm_tuple[g.single1]
    @inbounds permuted_idx = perm[phase_idxa]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[g.idx1*2-1] = bit1
    @inbounds phase[g.idx1*2] = bit2
    # second qubit permutation
    @inbounds perm = good_perm_tuple[g.single2]
    @inbounds permuted_idx = perm[phase_idxb]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[g.idx2*2-1] = bit1
    @inbounds phase[g.idx2*2] = bit2
    # bilateral
    @inbounds phase_idx = bit_to_int(phase[g.idx1*2-1],phase[g.idx1*2],phase[g.idx2*2-1], phase[g.idx2*2])
    @inbounds permuted_idx = cnot_perm[phase_idx]
    changed_phases = int_to_bit(permuted_idx, Val(4))
    @inbounds phase[g.idx1*2-1] = changed_phases[1]
    @inbounds phase[g.idx1*2] = changed_phases[2]
    @inbounds phase[g.idx2*2-1] = changed_phases[3]
    @inbounds phase[g.idx2*2] = changed_phases[4]
    return state
end

function QuantumClifford.apply!(state::BellState, g::GoodSingleQubitPerm) # TODO abstract away the permutation application as it is used by other gates too
    phase = state.phases
    # first qubit permutation
    @inbounds phase_idx = bit_to_int(phase[g.idx*2-1],phase[g.idx*2])
    if phase_idx==1
        return state
    end
    @inbounds perm = good_perm_tuple[g.single]
    @inbounds permuted_idx = perm[phase_idx]
    bit1, bit2 = int_to_bit(permuted_idx,Val(2))
    @inbounds phase[g.idx*2-1] = bit1
    @inbounds phase[g.idx*2] = bit2
    return state
end

##############################
# Noisy
##############################

"""A wrapper for `BellGate` that implements Pauli noise in addition to the gate."""
struct PauliNoiseBellGate{G} <: BellOp where {G<:BellOp} # TODO make it work with the QuantumClifford noise ops
    g::G
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.apply!(state::BellState, g::PauliNoiseBellGate)
    apply!(state, g.g)
    apply!(state, PauliNoiseOp(g.g.idx1,g.px,g.py,g.pz))
    apply!(state, PauliNoiseOp(g.g.idx2,g.px,g.py,g.pz))
    return state
end

"""`PauliNoiseOp(idx,px,py,pz)` causes qubit-pair `idx` to flip to one of the other 3 Bell states with probabilities `px`, `py`, `pz` respectively.

```jldoctest
julia> apply!(BellState([0,0]), PauliNoiseOp(1,1,0,0))
BellState(Bool[0, 1])

julia> apply!(BellState([0,0]), PauliNoiseOp(1,0,1,0))
BellState(Bool[1, 1])

julia> apply!(BellState([0,0]), PauliNoiseOp(1,0,0,1))
BellState(Bool[1, 0])
```
"""
struct PauliNoiseOp <: BellOp # TODO make it work with the QuantumClifford noise ops
    idx::Int
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.apply!(state::BellState, g::PauliNoiseOp)
    i = g.idx
    # TODO repetition with ...NoisyReset and PauliNoise...
    r = rand()
    if r<g.px
        apply!(state, BellPauliPermutation(2, i))
    elseif r<g.px+g.pz
        apply!(state, BellPauliPermutation(3, i))
    elseif r<g.px+g.pz+g.py
        apply!(state, BellPauliPermutation(4, i))
    end
    return state
end

"""Simulates twirled T1 noise"""
struct T1NoiseOp <: BellOp
    idx::Int
    λ₁::Float64
end

function QuantumClifford.apply!(state::BellState, g::T1NoiseOp)
    phase = state.phases
    input_state = bit_to_int(phase[g.idx*2-1],phase[g.idx*2]) # this is my initial state

    r = rand()
    λ₁ = g.λ₁

    output_state = if input_state==1
        if     r < 0.5*λ₁^2 - λ₁ + 1
            1
        elseif r < 0.5*λ₁^2 - λ₁ + 1  +  0.5*λ₁*(1-λ₁)
            2
        elseif r < 0.5*λ₁^2 - λ₁ + 1  +  0.5*λ₁*(1-λ₁)  +  0.5*λ₁^2
            3
        else # r < 1 = 0.5*λ₁^2 - λ₁ + 1  +  0.5*λ₁*(1-λ₁)  +  0.5*λ₁^2  +   0.5*λ₁*(1-λ₁)
            4
        end
    elseif input_state==2
        if     r < 0.5*λ₁
            1
        elseif r < 0.5*λ₁ + 1-λ₁
            2
        elseif r < 0.5*λ₁ + 1-λ₁  +  0.5*λ₁
            3
        else # r < 1 = 0.5*λ₁ + 1-λ₁  +  0.5*λ₁  + 0
            4
        end
    elseif input_state==3
        if     r < 0.25*λ₁*(λ₁+1)
            1
        elseif r < 0.25*λ₁*(λ₁+1)  +  0.5*λ₁*(1-λ₁)
            2
        elseif r < 0.25*λ₁*(λ₁+1)  +  0.5*λ₁*(1-λ₁)  +  0.25*λ₁*(λ₁+1)
            3
        else # r < 1 = 0.25*λ₁*(λ₁+1)  +  0.5*λ₁*(1-λ₁)  +  0.25*λ₁*(λ₁+1)  +  0.5*λ₁*(1-λ₁)
            4
        end
    else # input_state==4
        if     r < 0.5*λ₁
            1
        # elseif r < 0.5*λ₁ + 0             # output_state 2 is never reached
        #     2
        elseif r < 0.5*λ₁ + 0  +  0.5*λ₁
            3
        else # r < 1 = 0.5*λ₁ + 0  +  0.5*λ₁ + 1-λ₁
            4
        end
    end

    bit1, bit2 = int_to_bit(output_state, Val(2))
    @inbounds phase[g.idx*2-1] = bit1
    @inbounds phase[g.idx*2] = bit2
    return state
end

"""Simulates T2 noise"""
struct T2NoiseOp <: BellOp
    idx::Int
    λ₂::Float64
end

function QuantumClifford.apply!(state::BellState, g::T2NoiseOp)
    phase = state.phases
    input_state = bit_to_int(phase[g.idx*2-1],phase[g.idx*2]) # this is my initial state

    r = rand()
    λ₂ = g.λ₂

    output_state = if input_state==1
        if     r < 0.5*λ₂^2 - λ₂ + 1
            1
        else # r < 1 = 0.5*λ₂^2 - λ₂ + 1  +  0.5*λ₂*(2-λ₂)
            3
        end
    elseif input_state==2
        if     r < 0.5*λ₂^2 - λ₂ + 1
            2
        else # r < 1 = 0.5*λ₂^2 - λ₂ + 1  +  0.5*λ₂*(2-λ₂)
            4
        end
    elseif input_state==3
        if     r < 0.5*λ₂^2 - λ₂ + 1
            3
        else # r < 1 = 0.5*λ₂^2 - λ₂ + 1  +  0.5*λ₂*(2-λ₂)
            1
        end
    else # input_state==4
        if     r < 0.5*λ₂^2 - λ₂ + 1
            4
        else # r < 1 = 0.5*λ₂^2 - λ₂ + 1  +  0.5*λ₂*(2-λ₂)
            2
        end
    end

    bit1, bit2 = int_to_bit(output_state, Val(2))
    @inbounds phase[g.idx*2-1] = bit1
    @inbounds phase[g.idx*2] = bit2
    return state
end

"""A wrapper for [`BellMeasure`](@ref) that implements measurement noise."""
struct NoisyBellMeasure <: BellOp # TODO make it work with the QuantumClifford noise ops
    m::BellMeasure
    p::Float64
end

"""A wrapper for [`BellMeasure`](@ref) that implements measurement noise and Pauli noise after the reset."""
struct NoisyBellMeasureNoisyReset <: BellOp # TODO make it work with the QuantumClifford noise ops
    m::BellMeasure
    p::Float64
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.applywstatus!(state::BellState, op::NoisyBellMeasure)
    state, result = bellmeasure!(state, op.m)
    state, result⊻(rand()<op.p) ? continue_stat : failure_stat
end

function QuantumClifford.applywstatus!(state::BellState, op::NoisyBellMeasureNoisyReset)
    state, result = bellmeasure!(state, op.m)
    cont = result⊻(rand()<op.p)
    cont && apply!(state, PauliNoiseOp(op.m.sidx,op.px,op.py,op.pz))
    state, cont ? continue_stat : failure_stat
end

function NoisyBellSwap(idx1,idx2,px,py,pz)
    return PauliNoiseBellGate(BellSwap(idx1,idx2), px,py,pz)
end

##############################
# Random
##############################

"""
Random Bell diagonal state.
Input is the number of shared Bell pairs in the entanglement network.
"""
function Random.rand(::Type{BellState}, num_bell::Int)
    return BellState(BitArray(rand(Bool,num_bell*2)))
end

"""
Random [`BellGate`](@ref) on Bell pairs `i` and `j`.
"""
function Random.rand(::Type{BellGate}, i::Int,j::Int)
    return BellGate(rand(1:4),rand(1:4),rand(1:20),rand(1:6),rand(1:6),i,j)
end

"""
Random [`CNOTPerm`](@ref) on Bell pairs `i` and `j`.
"""
function Random.rand(::Type{CNOTPerm}, i::Int,j::Int)
    return CNOTPerm(rand(1:6),rand(1:6),i,j)
end

"""
Random [`GoodSingleQubitPerm`](@ref) on Bell pair `i`.
"""
function Random.rand(::Type{GoodSingleQubitPerm}, i::Int)
    return GoodSingleQubitPerm(rand(1:6),i)
end

"""
Random [`BellDoublePermutation`](@ref) on Bell pairs `i` and `j`.
"""
function Random.rand(::Type{BellDoublePermutation}, i::Int,j::Int)
    return BellDoublePermutation(rand(1:20),i,j)
end

"""
Random [`BellPauliPermutation`](@ref) on Bell pair `i`.
"""
function Random.rand(::Type{BellPauliPermutation}, i::Int)
    return BellPauliPermutation(rand(1:4),i)
end

"""
Random [`BellSinglePermutation`](@ref) on Bell pair `i`.
"""
function Random.rand(::Type{BellSinglePermutation}, i::Int)
    return BellSinglePermutation(rand(1:6),i)
end

"""
Initialize a random [`BellMeasure`](@ref) on Bell pair `i`.
"""
function Random.rand(::Type{BellMeasure}, i::Int)
    return BellMeasure(rand(1:3),i)
end

##############################
# Conversion from BP to QC
##############################

"""
Convert a Bell preserving gate from `BPGates` representation to a sequence of Clifford gate from `QuantumClifford`.

Translating [`BellMeasure`](@ref) is not precise,
because a real Bell measurement would destroy the Bell pair,
which can not be represented in `BPGates.jl`.
"""
function toQCcircuit end

function toQCcircuit(gate::BellSinglePermutation)
    return [
        SparseGate(one_perm_qc[gate.pidx], (gate.sidx*2-1,)),
        SparseGate(one_perm_qc[gate.pidx], (gate.sidx*2,))
    ]
end
function toQCcircuit(gate::BellDoublePermutation)
    return [
        SparseGate(double_perm_qc[gate.pidx], (gate.sidx1*2-1, gate.sidx2*2-1)),
        SparseGate(double_perm_qc[gate.pidx], (gate.sidx1*2, gate.sidx2*2))
    ]
end
function toQCcircuit(gate::BellPauliPermutation)
    return [pauli_perm_qc[gate.pidx](gate.sidx[1]*2-1)]
end

function toQCcircuit(gate::BellGate)
    return [
        pauli_perm_qc[gate.pauli1](gate.idx1*2-1),
        pauli_perm_qc[gate.pauli2](gate.idx2*2-1),
        SparseGate(double_perm_qc[gate.double], (gate.idx1*2-1, gate.idx2*2-1)),
        SparseGate(double_perm_qc[gate.double], (gate.idx1*2, gate.idx2*2)),
        SparseGate(one_perm_qc[gate.single1], (gate.idx1*2-1,)),
        SparseGate(one_perm_qc[gate.single1], (gate.idx1*2,)),
        SparseGate(one_perm_qc[gate.single2], (gate.idx2*2-1,)),
        SparseGate(one_perm_qc[gate.single2], (gate.idx2*2,)),
    ]
end

function toQCcircuit(gate::CNOTPerm)
    return [
        SparseGate(good_perm_qc[gate.single1][1], (gate.idx1*2-1,)),
        SparseGate(good_perm_qc[gate.single1][2], (gate.idx1*2,)),
        SparseGate(good_perm_qc[gate.single2][1], (gate.idx2*2-1,)),
        SparseGate(good_perm_qc[gate.single2][2], (gate.idx2*2,)),
        sCNOT(gate.idx1*2-1, gate.idx2*2-1),
        sCNOT(gate.idx1*2, gate.idx2*2)
    ]
end

function toQCcircuit(gate::GoodSingleQubitPerm)
    return [
        SparseGate(good_perm_qc[gate.single][1], (gate.idx*2-1,)),
        SparseGate(good_perm_qc[gate.single][2], (gate.idx*2,)),
    ]
end

function toQCcircuit(g::BellMeasure)
    meas = (sMX, sMY, sMZ)[g.midx]
    return [
        BellMeasurement([meas(g.sidx*2-1),meas(g.sidx*2)], g.midx==2)
        Reset(S"XX ZZ", [g.sidx*2-1,g.sidx*2])
    ]
end

function QuantumClifford.Stabilizer(state::BellState)
    res = bell(length(state.phases)÷2)
    tab(res).phases .= state.phases .* 0x2
    return res
end

function QuantumClifford.MixedDestabilizer(state::BellState)
    MixedDestabilizer(Stabilizer(state))
end

##############################
# Conversion from QC to BP
##############################

# TODO two pair gates
# TODO gates that do not preserve 0 and 00
# TODO general interface for all gates
function toBPpermutation1(circ) # this one works only on single pair mappings that preserve 0
    init_states = BellState.(Iterators.product(0:1,0:1))[:]
    init_stabs = Stabilizer.(init_states)
    final_stabs = [mctrajectory!(copy(s), circ)[1] for s in init_stabs]
    final_states = BellState.(final_stabs)
    init_ints = [bit_to_int(s.phases) for s in init_states]
    final_ints = [bit_to_int(s.phases) for s in final_states]
    i = findfirst(==(tuple(final_ints...)), good_perm_tuple)
    return GoodSingleQubitPerm(i,1)
end

function BellState(s::Stabilizer)
    r, c = size(s)::Tuple{Int, Int}
    r==c || throw(ArgumentError("Conversion to `BellState` failed. The stabilizer state has to be square in order to be convertible to a `BellState`."))
    s = canonicalize_rref!(copy(s))[1][end:-1:1]::Stabilizer
    bits = Bool[]
    for i in 1:r÷2
        j = (2i-1)
        s[j  ,j] == s[j  ,j+1] == (true, false) && all(==((false,false)), (s[j  ,k] for k in 1:c if k!=j && k!=j+1)) || throw(ArgumentError(lazy"Conversion to `BellState` failed. Row $(j  ) of the stabilizer state has to be of the form `..._XX_...` for it to be a valid `BellState`"))
        s[j+1,j] == s[j+1,j+1] == (false, true) && all(==((false,false)), (s[j+1,k] for k in 1:c if k!=j && k!=j+1)) || throw(ArgumentError(lazy"Conversion to `BellState` failed. Row $(j+1) row of the stabilizer state has to be of the form `..._ZZ_...` for it to be a valid `BellState`"))
        push!(bits, s.tab.phases[j]÷2)
        push!(bits, s.tab.phases[j+1]÷2)
    end
    return BellState(bits)
end

end # module
