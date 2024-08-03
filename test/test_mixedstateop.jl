using Test
using Logging

using BPGates: get_mixed_transition_probs, select_rand_state

using QuantumClifford: apply!

# TODO: have tolerance for numerical precision for noise calculations

# Set up the logger to enable debug level
function enable_debug_logging()
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

# Get the bell states used for testing
function get_bell_states()
    return (1 / sqrt(2)) * [
        [1.0; 0.0; 0.0; 1.0],
        [0.0; 1.0; 1.0; 0.0],
        [1.0; 0.0; 0.0; -1.0],
        [0.0; 1.0; -1.0; 0.0]
    ]
end

# Define kraus operators to test
function get_ampdamp_kraus_ops(λ::Float64)::Vector{Matrix{Complex}}
    # Define Kraus operators
    E1 = [1 0; 0 sqrt(1 - λ)]
    E2 = [0 sqrt(λ); 0 0]
    return [E1, E2]
end

# Returns an index selected from the indices of a list, weighted by the probability at that index
function get_weighted_index(prob_dist::Vector{Float64})
    if !isapprox(sum(prob_dist), 1.0, atol=1e-9)
        @error "get_weighted_index: prob_dist sums to $(sum(prob_dist)) and not 1"
    end
    rand_num = rand()
    cur_sum = 0.0
    for (idx, prob_val) in enumerate(prob_dist)
        if rand_num < cur_sum + prob_val
            return idx
        end
        cur_sum += prob_val
    end
end

# Converts a bit to an int from left to right (LSB is the rightmost bit)
function bit_to_int_testing(bits)
    reduce(⊻,(bit<<(length(bits) - index) for (index,bit) in enumerate(bits))) + 1 # +1 so that we use julia indexing convenctions
end

# Returns the tensor product of a set of Kraus operators
function get_tensored_kraus_ops(kraus_ops::Vector{Matrix{Complex}})
    tensored_kraus_ops = []
    for op1 in kraus_ops
        for op2 in kraus_ops
            push!(tensored_kraus_ops, kron(op1, op2))
        end
    end
    return tensored_kraus_ops
end

# Returns the diagonal of a matrix
function get_diagonal(matrix::Matrix{ComplexF64})::Vector{Float64}
    diagonal_length = min(size(matrix)[1], size(matrix)[2])
    if size(matrix)[1] != size(matrix)[2]
        @error("get_diagonal: matrix is not square, rows: $(size(matrix)[1]) != cols: $(size(matrix)[2])")
    end
    diag_elements = []
    for idx in 1:diagonal_length
        push!(diag_elements, matrix[idx, idx])
    end
    return diag_elements
end

# Returns the bell state corresponding to the input phases
function get_bell_state(phases::BitVector)
    state_idx = bit_to_int_testing(phases)
    bell_states = get_bell_states()
    return bell_states[state_idx]
end

# Numerically computes the transition probabilities of one state to another,
# given a mixed probability distribution of states and a set of kraus operators
function get_numerical_trans_probs(mixed_state_probs::Vector{Float64}, kraus_ops::Vector{Matrix{Complex}})
    bell_states = get_bell_states()
    ρ_bell = zeros(ComplexF64, length(bell_states), length(bell_states))
    for (bell_idx, bell_prob) in enumerate(mixed_state_probs)
        ρ_bell += bell_prob * bell_states[bell_idx] * bell_states[bell_idx]'
    end
    @debug("get_numerical_trans_probs: ρ_bell: $ρ_bell")
    tensored_ops = get_tensored_kraus_ops(kraus_ops)
    # Apply the Kraus operators to the Bell state
    ρ_new = zeros(ComplexF64, size(ρ_bell)[1], size(ρ_bell)[2])
    for K in tensored_ops
        ρ_new += K * ρ_bell * adjoint(K)
    end
    U_bell = hcat(get_bell_states()...)
    @debug("get_numerical_trans_probs: ρ_new: $ρ_new")
    @debug("get_numerical_trans_probs: U_bell: $U_bell")
    ρ_Bell_basis = adjoint(U_bell) * ρ_new * U_bell
    @debug("get_numerical_trans_probs: ρ_Bell_basis: $ρ_Bell_basis")
    @debug("get_numerical_trans_probs: size(ρ_Bell_basis): $(size(ρ_Bell_basis))")
    return get_diagonal(ρ_Bell_basis)
end

# Returns the simulated transition probabilities from the mixed state probabilities,
# the noise parameter, simulated for num_samples times
function get_simulated_trans_probs(mixed_state_probs::Vector{Float64}, noise_param::Float64, num_samples::Integer)
    sampled_state_count = fill(0, length(get_bell_states()))
    for _ in 1:num_samples
        state_idx = get_weighted_index(mixed_state_probs)
        transition_probs = get_mixed_transition_probs(noise_param, state_idx)
        new_state_idx = select_rand_state(transition_probs)
        # TODO: ensure that the index is inbounds for my list
        sampled_state_count[new_state_idx] += 1
    end
    return sampled_state_count / num_samples
end

# Testing functions

# Test our simulated noise for the pure states in the bell basis
function test_bell_basis_probs(kraus_ops_factory::Function)
    statedist_list = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    lambda = 1 - exp(-500e-9 / 260e-6)
    num_samples = 1000000
    kraus_ops = kraus_ops_factory(lambda)
    for state_dist in statedist_list
        @info "test_bell_basis_probs: state_dist: $state_dist"
        numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
        simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
        # @info("test_A_state_probs: numerical_trans_probs: $numerical_trans_probs")
        # @info("test_A_state_probs: simulated_trans_probs: $simulated_trans_probs")
        @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
    end
end

# Test our simulated noise when there should be no noise applied for states in the bell basis
function test_nonoise_bell_basis(kraus_ops_factory::Function)
    statedist_list = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    lambda = 0.0
    num_samples = 1000000
    kraus_ops = kraus_ops_factory(lambda)
    for state_dist in statedist_list
        @info "test_nonoise_bell_basis: state_dist: $state_dist"
        numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
        simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
        # @info("test_A_state_probs: numerical_trans_probs: $numerical_trans_probs")
        # @info("test_A_state_probs: simulated_trans_probs: $simulated_trans_probs")
        @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
    end
end

# Test our simulated noise with a noise factor close to 1 for states in the bell basis
function test_almost_all_noise(kraus_ops_factory::Function)
    statedist_list = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    lambda = 0.999999999
    num_samples = 1000000
    kraus_ops = kraus_ops_factory(lambda)
    for state_dist in statedist_list
        @info "test_almost_all_noise: state_dist: $state_dist"
        numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
        simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
        @info("test_almost_all_noise: numerical_trans_probs: $numerical_trans_probs")
        @info("test_almost_all_noise: simulated_trans_probs: $simulated_trans_probs")
        @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
    end
end

# Test random noise factors for states in the bell basis
function test_random_lambdas(kraus_ops_factory::Function)
    statedist_list = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    num_samples = 1000000
    num_random_lambdas = 5
    for _ in 1:num_random_lambdas
        lambda = rand()
        num_samples = 1000000
        kraus_ops = kraus_ops_factory(lambda)
        for state_dist in statedist_list
            @info "test_random_lambdas: state_dist: $state_dist"
            @info "test_random_lambdas: lambda: $lambda"
            numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
            simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
            @info("test_random_lambdas: numerical_trans_probs: $numerical_trans_probs")
            @info("test_random_lambdas: simulated_trans_probs: $simulated_trans_probs")
            @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
        end
    end
end

# Test all mixed states composed of two states in the Bell basis with equal probability
function test_equal_mixed_probs_2(kraus_ops_factory::Function)
    bell_states = get_bell_states()
    lambda = 1 - exp(-500e-9 / 260e-6)
    kraus_ops = kraus_ops_factory(lambda)
    num_samples = 1000000
    for first_phases_idx in 1:length(bell_states)
        for second_phases_idx in 1:length(bell_states)
            if first_phases_idx == second_phases_idx
                continue
            end
            state_dist = [0.0 for _ in 1:length(bell_states)]
            for basis_idx in 1:length(bell_states)
                if basis_idx == first_phases_idx || basis_idx == second_phases_idx
                    state_dist[basis_idx] = 0.5
                end
            end
            numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
            simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
            @info("test_equal_mixed_probs_2: numerical_trans_probs: $numerical_trans_probs")
            @info("test_equal_mixed_probs_2: simulated_trans_probs: $simulated_trans_probs")
            @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
        end
    end
end

# Test a mixed state composed of all 4 bell basis states with equal probability
function test_equal_mixed_probs_4(kraus_ops_factory::Function)
    state_dist = [0.25 for _ in 1:length(get_bell_states())]
    lambda = 1 - exp(-500e-9 / 260e-6)
    num_samples = 1000000
    kraus_ops = kraus_ops_factory(lambda)
    @info "test_equal_mixed_probs_4: state_dist: $state_dist"
    numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
    simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
    @info("test_equal_mixed_probs_4: numerical_trans_probs: $numerical_trans_probs")
    @info("test_equal_mixed_probs_4: simulated_trans_probs: $simulated_trans_probs")
    @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
end

# Test random mixed states in the bell basis
function test_random_mixed_states(kraus_ops_factory::Function)
    num_random_states = 20
    statedist_list = []
    for _ in 1:num_random_states
        state_dist = [rand() for _ in 1:length(get_bell_states())]
        state_dist = state_dist / sum(state_dist)
        push!(statedist_list, state_dist)
    end
    @debug "length(statedist_list): $(length(statedist_list))"
    lambda = 1 - exp(-500e-9 / 260e-6)
    num_samples = 1000000
    kraus_ops = kraus_ops_factory(lambda)
    for state_dist in statedist_list
        @info "test_random_mixed_states: state_dist: $state_dist"
        numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
        simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
        @info("test_random_mixed_states: numerical_trans_probs: $numerical_trans_probs")
        @info("test_random_mixed_states: simulated_trans_probs: $simulated_trans_probs")
        @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
    end
end

# Test random mixed states with random noise parameters in the Bell basis
function test_random_mixed_states_random_lambdas(kraus_ops_factory::Function)
    num_random_states = 20
    statedist_list = []
    for _ in 1:num_random_states
        state_dist = [rand() for _ in 1:length(get_bell_states())]
        state_dist = state_dist / sum(state_dist)
        push!(statedist_list, state_dist)
    end
    num_random_lambdas = 20
    num_samples = 1000000
    for _ in 1:num_random_lambdas
        lambda = rand()
        @info "test_random_mixed_states_random_lambdas: lambda: $lambda"
        kraus_ops = kraus_ops_factory(lambda)
        for state_dist in statedist_list
            @info "test_random_mixed_states_random_lambdas: state_dist: $state_dist"
            numerical_trans_probs = get_numerical_trans_probs(state_dist, kraus_ops)
            simulated_trans_probs = get_simulated_trans_probs(state_dist, lambda, num_samples)
            @info("test_random_mixed_states_random_lambdas: numerical_trans_probs: $numerical_trans_probs")
            @info("test_random_mixed_states_random_lambdas: simulated_trans_probs: $simulated_trans_probs")
            @test all(isapprox.(numerical_trans_probs, simulated_trans_probs, atol=(2 / sqrt(num_samples))))
        end
    end
end

# Run all tests
@testset "BPGates.jl, get_mixed_transition_probs tests" begin
    enable_debug_logging()
    test_bell_basis_probs(get_ampdamp_kraus_ops)
    test_nonoise_bell_basis(get_ampdamp_kraus_ops)
    test_almost_all_noise(get_ampdamp_kraus_ops)
    test_random_lambdas(get_ampdamp_kraus_ops)
    test_equal_mixed_probs_2(get_ampdamp_kraus_ops)
    test_equal_mixed_probs_4(get_ampdamp_kraus_ops)
    test_random_mixed_states(get_ampdamp_kraus_ops)
    test_random_mixed_states_random_lambdas(get_ampdamp_kraus_ops)
end
