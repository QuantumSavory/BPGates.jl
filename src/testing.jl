using Revise
using QuantumClifford
using BPGates
circuit1 = [
    BellSinglePermutation(1, 1),
    BellSinglePermutation(2, 2),
    BellSinglePermutation(3, 1),   # pair 1 used again -> pair 2 idle in between
]

model1 = BPCircuitNoise(
    single_pair = PauliNoiseData(0.01, 0.01, 0.01),
    idle_noise  = PauliNoiseData(0.001, 0.001, 0.001),
)
 
result1 = noisify(circuit1, model1)
for op in result1
    println(typeof(op), "  ", op)
end

