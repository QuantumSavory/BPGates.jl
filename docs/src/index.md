# BPGates.jl

```@meta
DocTestSetup = quote
    using QuantumClifford, BPGates
end
```

Faster Bell-preserving gates for Clifford circuit simulators like [QuantumClifford.jl](https://github.com/Krastanov/QuantumClifford.jl).

```
using BPGates
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using Random

npairs = 10
state = BellState(npairs)
random_gates = [rand(BellGate,randperm(npairs)[1:2]...) for _ in 1:10]
px = py = pz = 0.01
noisy_random_gates = [PauliNoiseBellGate(g,px,py,pz) for g in random_gates]
random_measurements = [rand(BellMeasure,rand(1:npairs)) for _ in 1:10]
p = 0.01
px0 = py0 = pz0 = 0.03
noisy_random_measurements = [NoisyBellMeasureNoisyReset(m,p,px,py,pz) for m in random_measurements]
all_ops = vcat(noisy_random_gates,noisy_random_measurements)
random_circuit = all_ops[randperm(npairs)]
mctrajectory!(copy(state), random_circuit)
# TODO initialization network noise
# TODO measurement noise
# TODO more restricted types of good gates
```