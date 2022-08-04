using QuantumClifford.Experimental.NoisyCircuits
using BPGates: convert2QC
using Random

function test_viaquantumclifford()
    @testset "Comparing to QuantumClifford" begin
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
end

test_viaquantumclifford()
