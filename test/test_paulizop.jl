using Test
using Logging

using BPGates: get_mixed_transition_probs, select_rand_state

using .TestKrausNoiseHelpers

using QuantumClifford: apply!

# TODO: have tolerance for numerical precision for noise calculations

# Set up the logger to enable debug level
function enable_debug_logging()
    global_logger(ConsoleLogger(stderr, Logging.Info))
end

# Define kraus operators to test
function get_phasedamp_kraus_ops(λ::Float64)::Vector{Matrix{Complex}}
    # Define Kraus operators
    E1 = sqrt(1 - (λ / 2)) * [1 0; 0 1]
    E2 = sqrt(λ / 2) * [1 0; 0 -1]
    return [E1, E2]
end

function do_z_flip(state_idx::Int)
    if state_idx == 1
        return 3
    elseif state_idx == 2
        return 4
    elseif state_idx == 3
        return 1
    elseif state_idx == 4
        return 2
    else
        @error "do_z_flip: state_idx is not between 1 and 4, state_idx: $state_idx"
    end
end

# Returns the simulated transition probabilities from the mixed state probabilities,
# the noise parameter, simulated for num_samples times
function get_simulated_trans_probs(mixed_state_probs::Vector{Float64}, noise_param::Float64, num_samples::Integer)
    sampled_state_count = fill(0, length(TestKrausNoiseHelpers.get_bell_states()))
    for _ in 1:num_samples
        state_idx = TestKrausNoiseHelpers.get_weighted_index(mixed_state_probs)
        # transition_probs = TestKrausNoiseHelpers.get_mixed_transition_probs(noise_param, state_idx)
        # new_state_idx = TestKrausNoiseHelpers.select_rand_state(transition_probs)
        r = rand()
        if r < 2 * (noise_param / 2) * (1 - (noise_param / 2))
            new_state_idx = do_z_flip(state_idx)
        else
            new_state_idx = state_idx
        end
        # TODO: ensure that the index is inbounds for my list
        sampled_state_count[new_state_idx] += 1
    end
    return sampled_state_count / num_samples
end

# Run all tests
@testset "BPGates.jl, get_mixed_transition_probs tests" begin
    enable_debug_logging()
    test_bell_basis_probs(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_nonoise_bell_basis(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_almost_all_noise(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_random_lambdas(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_equal_mixed_probs_2(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_equal_mixed_probs_4(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_random_mixed_states(get_phasedamp_kraus_ops, get_simulated_trans_probs)
    test_random_mixed_states_random_lambdas(get_phasedamp_kraus_ops, get_simulated_trans_probs)
end
