@testitem "noisify" begin
    using QuantumClifford: affectedqubits, noisify, AbstractNoise, PauliNoise, CircuitNoise
    using BPGates: T1Noise, T2Noise, T1T2Noise, BellMeasure, MeasurementFlipNoise

    @testset "CNOTPerm gate noise applies only to gate-affected slots" begin
        g = CNOTPerm(1, 1, 2, 5)
        q1, q2 = affectedqubits(g)

        result = noisify(g, CircuitNoise(two_qubit = PauliNoise(0.1, 0.1, 0.1)))

        @test result[1] === g
        noise_ops = result[2:end]
        @test length(noise_ops) == 2
        @test Set(op.idx for op in noise_ops) == Set([q1, q2])
        @test all(op -> op isa PauliNoiseOp, noise_ops)
    end

    @testset "idle slot tracking correct across mixed op types" begin
        circuit = [BellSwap(1, 2), CNOTPerm(2, 2, 1, 3), BellMeasure(1, 2)]
        result = noisify(circuit, CircuitNoise(idle_noise = T1Noise(0.1)))

        noise_ops = filter(op -> op isa T1NoiseOp, result)
        @test length(noise_ops) == 1
        @test noise_ops[1].idx == 3
    end

    @testset "T1T2Noise idle noise inserts both T1 and T2 ops per idle slot" begin
        # Unlike T1Noise/T2Noise alone, T1T2Noise represents both mechanisms
        # acting during the same idle window: append_idle_noise! pushes a
        # T1NoiseOp immediately followed by a T2NoiseOp for each idle gap on
        # a given qubit -- NOT a single combined op type.
        circuit = [BellSwap(1, 2), CNOTPerm(2, 2, 1, 3), BellMeasure(1, 2)]
        result = noisify(circuit, CircuitNoise(idle_noise = T1T2Noise(0.1, 0.2)))

        t1_ops = filter(op -> op isa T1NoiseOp, result)
        t2_ops = filter(op -> op isa T2NoiseOp, result)

        @test length(t1_ops) == 1
        @test length(t2_ops) == 1
        @test t1_ops[1].idx == 3
        @test t2_ops[1].idx == 3
        @test t1_ops[1].λ₁ == 0.1
        @test t2_ops[1].λ₂ == 0.2

        # confirm ordering: T1NoiseOp immediately precedes its paired T2NoiseOp
        t1_pos = findfirst(op -> op isa T1NoiseOp, result)
        t2_pos = findfirst(op -> op isa T2NoiseOp, result)
        @test t2_pos == t1_pos + 1
    end

    @testset "measurement noise applies flip probability at correct slot" begin
        m = BellMeasure(1, 2)   # midx=1 (basis), sidx=2 (the actual slot)

        noisy_result = noisify(m, CircuitNoise(measurement = MeasurementFlipNoise(0.05)))
        @test noisy_result == [NoisyBellMeasureNoisyReset(m, 0.05, 0.0, 0.0, 0.0)]

        clean_result = noisify(m, CircuitNoise())
        @test clean_result == [m]
    end

    @testset "already-noisy ops pass through unchanged, even with new noise configured" begin
        # Per the noisify issue spec: "NoiseOp and NoisyGate should not
        # accidentally get double-noisified unless explicitly requested."
        # This is a uniform rule across every already-noisy op type --
        # a second CircuitNoise pass never modifies an op that is already
        # noisy, regardless of type. No compose-on-second-pass exceptions.
        noise = CircuitNoise(
            single_qubit = T1Noise(0.1),
            two_qubit    = PauliNoise(0.2, 0.2, 0.2),
            measurement  = MeasurementFlipNoise(0.3),
        )

        t1   = T1NoiseOp(2, 0.1)
        t2   = T2NoiseOp(2, 0.1)
        pn   = PauliNoiseOp(1, 0.1, 0.1, 0.1)
        nbm  = NoisyBellMeasure(BellMeasure(1, 2), 0.1)
        nbmr = NoisyBellMeasureNoisyReset(BellMeasure(1, 2), 0.1, 0.0, 0.0, 0.0)

        @test noisify(t1,   noise) == [t1]
        @test noisify(t2,   noise) == [t2]
        @test noisify(pn,   noise) == [pn]
        @test noisify(nbm,  noise) == [nbm]
        @test noisify(nbmr, noise) == [nbmr]
    end

    @testset "no noise configured leaves CNOTPerm/BellMeasure circuit untouched" begin
        circuit = [
            BellGate(1, 2, 1, 2, 1, 2, 1),
            CNOTPerm(1, 1, 1, 3),
            BellMeasure(1, 2),
        ]
        result = noisify(circuit, CircuitNoise())
        @test result == circuit
        @test all(op -> !(op isa Union{T1NoiseOp,T2NoiseOp,PauliNoiseOp}), result)
        @test !any(op -> op isa NoisyBellMeasureNoisyReset, result)
    end
end