using Test
using Logging

using BPGates

using QuantumClifford

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
    
    QuantumClifford.apply!(test_state, test_guaranteed_pauliz)

    @test test_state.phases[2*changed_bp_idx - 1] == 0
    @test test_state.phases[2*changed_bp_idx] == 1
end

function test_pauliz_application_guaranteed_none()
    n_bellpairs = 3
    changed_bp_idx = 2
    test_state = BellState(n_bellpairs)
    test_guaranteed_pauliz = PauliZOp(changed_bp_idx, 0)
    
    # TODO: do I have to import QuantumClifford to use it when testing?
    QuantumClifford.apply!(test_state, test_guaranteed_pauliz)

    @test test_state.phases[2*changed_bp_idx - 1] == 0
    @test test_state.phases[2*changed_bp_idx] == 0
end

function test_mixed_state_op_constructor()
    mixed_state_op = MixedStateOp(1, 0.5)
    
    @test mixed_state_op isa MixedStateOp

    @test mixed_state_op.idx == 1

    @test mixed_state_op.lambda == 0.5
end

function test_apply_mixed_state_op()
    n_bellpairs = 1
    changed_bp_idx = 1
    test_state = BellState(n_bellpairs)
    mixed_state_op = MixedStateOp(changed_bp_idx, 0.0)

    QuantumClifford.apply!(test_state, mixed_state_op)

    @test test_state.phases[2 * changed_bp_idx - 1] == 0
    @test test_state.phases[2 * changed_bp_idx] == 0
end

function test_apply_mixed_state_op_diff_bellstate()
    n_bellpairs = 1
    changed_bp_idx = 1
    test_state = BellState((0, 1))
    mixed_state_op = MixedStateOp(changed_bp_idx, 0.0)
    
    QuantumClifford.apply!(test_state, mixed_state_op)
    
    @test test_state.phases[2 * changed_bp_idx - 1] == 0
    @test test_state.phases[2 * changed_bp_idx] == 1
end

function test_apply_both_memory_errors()
    n_bellpairs = 1
    changed_bp_idx = 1
    
    test_state = BellState(n_bellpairs)

    noise_ops = Vector{AbstractNoiseBellOp}()
    
    push!(noise_ops, PauliZOp(changed_bp_idx, 0))
    push!(noise_ops, MixedStateOp(changed_bp_idx, 0.0))

    for noise_op in noise_ops
        QuantumClifford.apply!(test_state, noise_op)
    end

    @test test_state.phases[2 * changed_bp_idx - 1] == 0
    @test test_state.phases[2 * changed_bp_idx] == 0
end

@testset "BPGates.jl, PauliZOp tests" begin
    test_pauliz_constructor()
    test_pauliz_application_guaranteed()
    test_pauliz_application_guaranteed_none()
end

@testset "BPGates.jl, MixedStateOp tests" begin
    test_mixed_state_op_constructor()
    test_apply_mixed_state_op()
    test_apply_mixed_state_op_diff_bellstate()
    test_apply_both_memory_errors()
end
