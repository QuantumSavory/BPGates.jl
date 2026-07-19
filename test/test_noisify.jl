@testitem "noisify" begin
    using QuantumClifford: affectedqubits, noisify, AbstractNoise, PauliNoise, CircuitNoise
    using BPGates: T1Noise, T2Noise, BellMeasure, MeasurementFlipNoise

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

    @testset "measurement noise applies flip probability at correct slot" begin
        m = BellMeasure(1, 2)   # midx=1 (basis), sidx=2 (the actual slot)

        noisy_result = noisify(m, CircuitNoise(measurement = MeasurementFlipNoise(0.05)))
        @test noisy_result == [NoisyBellMeasureNoisyReset(m, 0.05, 0.0, 0.0, 0.0)]

        clean_result = noisify(m, CircuitNoise())
        @test clean_result == [m]
    end

    @testset "measurement noise composes while preserving type" begin
        noise = CircuitNoise(measurement = MeasurementFlipNoise(0.2))
        expected_p = 0.1 + 0.2 - 2 * 0.1 * 0.2

        nbm = NoisyBellMeasure(BellMeasure(1, 2), 0.1)
        result_nbm = noisify(nbm, noise)

        @test result_nbm == [
            NoisyBellMeasure(nbm.m, expected_p)
        ]
        @test only(result_nbm) isa NoisyBellMeasure

        nbmr = NoisyBellMeasureNoisyReset(
            BellMeasure(1, 2),
            0.1,
            0.0,
            0.0,
            0.0,
        )

        result_nbmr = noisify(nbmr, noise)

        @test result_nbmr == [
            NoisyBellMeasureNoisyReset(
                nbmr.m,
                expected_p,
                nbmr.px,
                nbmr.py,
                nbmr.pz,
            )
        ]
        @test only(result_nbmr) isa NoisyBellMeasureNoisyReset

        @test noisify(nbm, CircuitNoise()) == [nbm]
        @test noisify(nbmr, CircuitNoise()) == [nbmr]
    end

    @testset "already-noisy gate ops pass through unchanged, even with new noise configured" begin
        # NOTE: PauliNoiseOp/T1NoiseOp/T2NoiseOp do not currently compose additional
        # noise the way NoisyBellMeasureNoisyReset does — a second CircuitNoise pass
        # is silently a no-op for these types. This is a known, intentional gap for
        # now (see PR discussion), not yet fixed to append a second noise op.
        noise = CircuitNoise(single_qubit = T1Noise(0.1), two_qubit = PauliNoise(0.2, 0.2, 0.2))

        t1 = T1NoiseOp(2, 0.1)
        @test noisify(t1, noise) == [t1]
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