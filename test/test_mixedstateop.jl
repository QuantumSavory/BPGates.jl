using .TestKrausNoiseHelpers
using Logging

# TODO: have tolerance for numerical precision for noise calculations

# Set up the logger to enable debug level
function enable_debug_logging()
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

# Define kraus operators to test
function get_ampdamp_kraus_ops(λ::Float64)::Vector{Matrix{Complex}}
    # Define Kraus operators
    E1 = [1 0; 0 sqrt(1 - λ)]
    E2 = [0 sqrt(λ); 0 0]
    return [E1, E2]
end

# Returns the simulated transition probabilities from the mixed state probabilities,
# the noise parameter, simulated for num_samples times
function get_simulated_trans_probs(mixed_state_probs::Vector{Float64}, noise_param::Float64, num_samples::Integer)
    sampled_state_count = fill(0, length(TestKrausNoiseHelpers.get_bell_states()))
    for _ in 1:num_samples
        state_idx = TestKrausNoiseHelpers.get_weighted_index(mixed_state_probs)
        transition_probs = TestKrausNoiseHelpers.get_mixed_transition_probs(noise_param, state_idx)
        new_state_idx = TestKrausNoiseHelpers.select_rand_state(transition_probs)
        # TODO: ensure that the index is inbounds for my list
        sampled_state_count[new_state_idx] += 1
    end
    return sampled_state_count / num_samples
end

# Run all tests
@testset "BPGates.jl, get_mixed_transition_probs tests" begin
    enable_debug_logging()
    test_bell_basis_probs(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_nonoise_bell_basis(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_almost_all_noise(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_random_lambdas(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_equal_mixed_probs_2(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_equal_mixed_probs_4(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_random_mixed_states(get_ampdamp_kraus_ops, get_simulated_trans_probs)
    test_random_mixed_states_random_lambdas(get_ampdamp_kraus_ops, get_simulated_trans_probs)
end
