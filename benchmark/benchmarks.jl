using BenchmarkTools
using StableRNGs
using Random
using BPGates
using QuantumClifford

const SUITE = BenchmarkGroup()

"""Return two distinct random indices in 1:n."""
function twoidx(rng, n)
    i = rand(rng, 1:n)
    j = rand(rng, 1:n-1)
    j >= i && (j += 1)
    return i, j
end

##############################
# apply! benchmarks
##############################

SUITE["apply"] = BenchmarkGroup(["apply"])

# Single-pair gates
SUITE["apply"]["BellPauliPermutation"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["BellPauliPermutation"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [BellPauliPermutation(rand(rng, 1:4), rand(rng, 1:$n)) for _ in 1:$n]
    ) evals=1
end

SUITE["apply"]["BellSinglePermutation"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["BellSinglePermutation"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [BellSinglePermutation(rand(rng, 1:6), rand(rng, 1:$n)) for _ in 1:$n]
    ) evals=1
end

# Two-pair gates
SUITE["apply"]["BellDoublePermutation"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["BellDoublePermutation"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [BellDoublePermutation(rand(rng, 1:20), twoidx(rng, $n)...) for _ in 1:$n]
    ) evals=1
end

SUITE["apply"]["BellGate"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["BellGate"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [rand(rng, BellGate, twoidx(rng, $n)...) for _ in 1:$n]
    ) evals=1
end

SUITE["apply"]["CNOTPerm"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["CNOTPerm"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [CNOTPerm(rand(rng, 1:6), rand(rng, 1:6), twoidx(rng, $n)...) for _ in 1:$n]
    ) evals=1
end

SUITE["apply"]["GoodSingleQubitPerm"] = BenchmarkGroup()
for n in [10, 100, 1000]
    SUITE["apply"]["GoodSingleQubitPerm"]["n$n"] = @benchmarkable begin
        for op in circuit
            apply!(state, op)
        end
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = [GoodSingleQubitPerm(rand(rng, 1:6), rand(rng, 1:$n)) for _ in 1:$n]
    ) evals=1
end

##############################
# mctrajectory! benchmarks
##############################

SUITE["mctrajectory"] = BenchmarkGroup(["mctrajectory"])

# Mixed circuit with various gate types (no measurements, to avoid early termination)
SUITE["mctrajectory"]["mixed_gates"] = BenchmarkGroup()
for n in [10, 100, 1000]
    num_gates = min(n, 100)
    SUITE["mctrajectory"]["mixed_gates"]["n$n"] = @benchmarkable begin
        mctrajectory!(state, circuit)
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = vcat(
            [rand(rng, BellGate, twoidx(rng, $n)...) for _ in 1:$num_gates],
            [BellPauliPermutation(rand(rng, 1:4), rand(rng, 1:$n)) for _ in 1:$num_gates],
            [BellSinglePermutation(rand(rng, 1:6), rand(rng, 1:$n)) for _ in 1:$num_gates],
        )
    ) evals=1
end

# Circuit with measurements (may terminate early on failure)
SUITE["mctrajectory"]["with_measurements"] = BenchmarkGroup()
for n in [10, 100]
    num_gates = n
    SUITE["mctrajectory"]["with_measurements"]["n$n"] = @benchmarkable begin
        mctrajectory!(state, circuit)
    end setup=(
        rng = StableRNG(42);
        state = BellState($n);
        circuit = vcat(
            [rand(rng, BellGate, twoidx(rng, $n)...) for _ in 1:$num_gates],
            [BellMeasure(rand(rng, 1:3), rand(rng, 1:$n)) for _ in 1:$num_gates],
        )
    ) evals=1
end

# Noisy circuit
SUITE["mctrajectory"]["noisy"] = BenchmarkGroup()
for n in [10, 100, 1000]
    num_gates = min(n, 100)
    SUITE["mctrajectory"]["noisy"]["n$n"] = @benchmarkable begin
        mctrajectory!(state, circuit)
    end setup=(
        rng = StableRNG(42);
        state = rand(rng, BellState, $n);
        circuit = vcat(
            [PauliNoiseBellGate(rand(rng, BellGate, twoidx(rng, $n)...), 0.01, 0.01, 0.01) for _ in 1:$num_gates],
            [PauliNoiseOp(rand(rng, 1:$n), 0.01, 0.01, 0.01) for _ in 1:$num_gates],
        )
    ) evals=1
end
