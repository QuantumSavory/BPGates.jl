using Test
using Logging

using BPGates

function test_pauliz_constructor()

    test_pauliz = PauliZOp(1, 0.5)
    
    @test test_pauliz isa PauliZOp

    @test test_pauliz.idx == 1

    @test test_pauliz.pz == 0.5
end

function test_pauliz_application_guaranteed()
    n_bellpairs = 3
    changed_bp_idx = 2
    test_state = BellState(n_bellpairs)
    test_guaranteed_pauliz = PauliZOp(changed_bp_idx, 1)
    
    apply!(test_state, test_guaranteed_pauliz)

    @test test_state[2*changed_bp_idx - 1] == 0
    @test test_state[2*changed_bp_idx] == 1
end

function test_pauliz_application_guaranteed_none()
    n_bellpairs = 3
    changed_bp_idx = 2
    test_state = BellState(n_bellpairs)
    test_guaranteed_pauliz = PauliZOp(changed_bp_idx, 0)
    
    # TODO: do I have to import QuantumClifford to use it when testing?
    QuantumClifford.apply!(test_state, test_guaranteed_pauliz)

    @test test_state[2*changed_bp_idx - 1] == 0
    @test test_state[2*changed_bp_idx] == 0
end

@testset "BPGates.jl, PauliZOp tests" begin
    test_pauliz_constructor()
    test_pauliz_application_guaranteed()
    test_pauliz_application_guaranteed_none()
end
