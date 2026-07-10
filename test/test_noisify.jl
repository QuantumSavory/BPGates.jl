@testitem "noisify" begin
    using QuantumClifford: affectedqubits, noisify, skip_idling_noise, append_idle_noise!,
                            AbstractNoise, PauliNoise
    using BPGates: T1Noise, T2Noise, BPCircuitNoise, BellMeasure
    # ──────────────────────────────────────────────────────────────────────────
    # 1. Gate noise on CNOTPerm applies to both affected slots
    # ──────────────────────────────────────────────────────────────────────────
    @testset "CNOTPerm gate noise applies only to gate-affected slots" begin
        g = CNOTPerm(1, 1, 2, 5)   # idx1=2, idx2=5
        q1, q2 = affectedqubits(g)

        result = noisify(g, BPCircuitNoise(gate_noise = T1Noise(0.1)))

        @test result[1] === g
        noise_ops = result[2:end]
        @test length(noise_ops) == 2
        @test Set(op.idx for op in noise_ops) == Set([q1, q2])
        @test all(op -> op isa T1NoiseOp && op.λ₁ == 0.1, noise_ops)
    end

    @testset "idle slot tracking correct across mixed op types" begin

        circuit = [BellSwap(1, 2), CNOTPerm(2, 2, 1, 3), BellMeasure(1, 2)]
        result = noisify(circuit, BPCircuitNoise(idle_noise = T1Noise(0.1)))

        noise_ops = filter(op -> op isa T1NoiseOp, result)
        @test length(noise_ops) == 1
        @test noise_ops[1].idx == 3
    end
    @testset "measurement noise applies flip probability at correct slot" begin
        m = BellMeasure(1, 2)   # midx=1 (basis), sidx=2 (the actual slot)

        noisy_result = noisify(m, BPCircuitNoise(measurement = 0.05))
        @test noisy_result == [NoisyBellMeasureNoisyReset(m, 0.05, 0.0, 0.0, 0.0)]

        clean_result = noisify(m, BPCircuitNoise())
        @test clean_result == [m]
    end
    @testset "already-noisy ops pass through unchanged" begin
        noise = BPCircuitNoise(gate_noise = T1Noise(0.1), measurement = 0.2)

        t1   = T1NoiseOp(2, 0.1)
        nbmr = NoisyBellMeasureNoisyReset(BellMeasure(1, 2), 0.2, 0.0, 0.0, 0.0)

        @test noisify(t1,   noise) == [t1]
        @test noisify(nbmr, noise) == [nbmr]
    end

    @testset "no noise configured leaves CNOTPerm/BellMeasure circuit untouched" begin
        circuit = [
            BellGate(1, 2, 1, 2, 1, 2, 1),
            CNOTPerm(1, 1, 1, 3),
            BellMeasure(1, 2),
        ]
        result = noisify(circuit, BPCircuitNoise())
        print(result)
        @test result == circuit
        @test all(op -> !(op isa Union{T1NoiseOp,T2NoiseOp,PauliNoiseOp}), result)
        @test !any(op -> op isa NoisyBellMeasureNoisyReset, result)
    end
end