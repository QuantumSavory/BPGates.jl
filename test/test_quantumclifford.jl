using QuantumClifford.Experimental.NoisyCircuits
using BPGates: convert2QC
using Random

function test_viaquantumclifford()
    @testset "BPGate comparison" begin
        for num_bell in test_sizes[2:end]
            num_gates = 10
            state = BellState(num_bell)
            circuit = [rand(BellGate,randperm(num_bell)[1:2]...) for _ in 1:num_gates]

            endstate, status = mctrajectory!(copy(state), circuit)

            stabstate = Stabilizer(copy(state))
            stabcircuit = vcat([convert2QC(g) for g in circuit]...) # TODO make the conversion neater
            for (g,idx) in stabcircuit
                apply!(stabstate,g,collect(idx))
            end
            @test canonicalize!(stabstate) == canonicalize!(Stabilizer(endstate))
        end
    end
    @testset "CNOTPerm comparison" begin
        for num_bell in test_sizes[2:end]
            num_gates = 10
            state = BellState(num_bell)
            circuit = [CNOTPerm(1,1,randperm(num_bell)[1:2]...) for _ in 1:num_gates]

            endstate, status = mctrajectory!(copy(state), circuit)

            stabstate = Stabilizer(copy(state))
            stabcircuit = vcat([[sCNOT(g.idx1*2-1, g.idx2*2-1), sCNOT(g.idx1*2, g.idx2*2)] for g in circuit]...) # TODO make the conversion neater
            for g in stabcircuit
                apply!(stabstate,g)
            end
            @test canonicalize!(stabstate) == canonicalize!(Stabilizer(endstate))
        end
    end
end

test_viaquantumclifford()
