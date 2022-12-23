# [Noisy Gates](@id noisy-gates)

```@meta
DocTestSetup = quote
    using QuantumClifford, BPGates
end
```

## Pauli Noise

[`PauliNoiseOp`](@ref) performs Pauli noise of given magnitude (chaces for X, Y, and Z errors). It is used as a primitive inside of noisy operations.

```jldoctest
julia> apply!(BellState([0,0]), PauliNoiseOp(1,1,0,0))
BellState(Bool[0, 1])

julia> apply!(BellState([0,0]), PauliNoiseOp(1,0,1,0))
BellState(Bool[1, 1])

julia> apply!(BellState([0,0]), PauliNoiseOp(1,0,0,1))
BellState(Bool[1, 0])
```

## Noisy gate wrapper

[`PauliNoiseBellGate`](@ref) can be used to wrap a normal gate with a noise process.

## Noisy measurements

[`NoisyBellMeasureNoisyReset`](@ref) performs a `BellMeasurement` with a chance `p` to report the opposite result for a coincidence measurement, and chances `px`, `py`, and `pz` to flip the new Bell state to (respectively) one of the other 3 Bell states.

It implements `applywstatus!` which enables its use with `mctrajectories` from `QuantumClifford`.