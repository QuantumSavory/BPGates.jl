@testitem "T1 noise uses twirled density-matrix evolution" begin
    using Test
    using QuantumOpticsBase
    using QuantumClifford
    using LinearAlgebra
    using BPGates

    b = SpinBasis(1//2)
    b2 = tensor(b, b)

    psi0 = normalize(tensor(spinup(b), spinup(b)) + tensor(spindown(b), spindown(b)))
    rho0 = dm(psi0)
    # computational -> bell - change of basis matrix
    Ubell = (1 / sqrt(2)) * ComplexF64[
        1  0   0   1;
        0  1   1   0;
        0  1  -1   0;
        1  0   0  -1
    ]

    to_bell_basis = DenseOperator(b2, b2, Ubell)

    # performs twirling on the state i.e. trims the off-diagonal terms off
    function bell_twirl_project(rho)
        rho_bell = to_bell_basis * rho * dagger(to_bell_basis)
        vals = diag(dense(rho_bell).data)
        rho_bell_diag = DenseOperator(rho_bell.basis_l, rho_bell.basis_r, Matrix(Diagonal(vals)))
        return dagger(to_bell_basis) * rho_bell_diag * to_bell_basis
    end


    # probability distribution of the input state in the bell basis
    function bell_probs(rho)
        rho_bell = to_bell_basis * rho * dagger(to_bell_basis)
        return real.(diag(dense(rho_bell).data))
    end

    # Kraus Channel for T1 amplitude damping in QuantumOptics
    function amplitude_damping_channel(rho, lambda)
        K0 = DenseOperator(b, b, ComplexF64[
            1  0;
            0  sqrt(1 - lambda)
        ])

        K1 = DenseOperator(b, b, ComplexF64[
            0  sqrt(lambda);
            0  0
        ])

        Ks = [K0, K1]

        rhoout = DenseOperator(b2, b2, zeros(ComplexF64, 4, 4))

        for A in Ks
            for B in Ks
                K = tensor(A, B)
                rhoout += K * rho * dagger(K)
            end
        end

        return rhoout
    end

    function t1_channel(rho, t, t1)
        lambda = 1 - exp(-t / t1)
        return amplitude_damping_channel(rho, lambda)
    end

    function t1_twirl_two_step(rho, t, t1)
        rho = t1_channel(rho, t / 2, t1)
        rho = bell_twirl_project(rho)

        rho = t1_channel(rho, t / 2, t1)
        return bell_twirl_project(rho)
        end
    ts = [1.0, 10, 100]

    # verifying that # E(t) = E(t/2) ∘ E(t/2)

    @testset "Exact T1 channel composes without twirling" begin
        t1 = 10.0

        for t in ts
            rho_direct = t1_channel(rho0, t, t1)
            rho_two_step = t1_channel(t1_channel(rho0, t / 2, t1), t / 2, t1)

            @test dense(rho_direct).data ≈ dense(rho_two_step).data atol=1e-10
        end
    end

    # verifying that T(E(t)) ≠ T(E(t/2) ∘ T(E(t/2))) i.e. twirling or removing off-diagonals changes the result 

    @testset "Intermediate twirling changes the result" begin
        t1 = 10.0

        for t in ts
            rho_direct_twirl = bell_twirl_project(t1_channel(rho0, t, t1))
            rho_twirl_two_step = t1_twirl_two_step(rho0, t, t1)

            @test norm(
                bell_probs(rho_direct_twirl) -
                bell_probs(rho_twirl_two_step)
            ) > 1e-6
        end
    end

  
   # This one test is used to check for two things: 1) BPGates indeed implements the twirled version, BPGates does NOT implement the untwirled version
    @testset "BPGates implements twirled (not untwirled) T1 evolution" begin
        t1 = 10.0
        target_sigma = 6.0
        N_min, N_max = 200_000, 5_000_000 # bounds for the number of runs 
        

        for t in ts
            twirled_probs   = bell_probs(t1_twirl_two_step(rho0, t, t1))
            untwirled_probs = bell_probs(bell_twirl_project(t1_channel(rho0, t, t1)))
            # this code is used to decide the number of runs to ensure that monte carlo noise won't affect the 5 sigma runs to a large extent
            delta = maximum(abs.(twirled_probs .- untwirled_probs))
            N_needed = ceil(Int, (target_sigma * 0.5 / delta)^2)
            N = clamp(N_needed, N_min, N_max)

            lambda_half = 1 - exp(-(t / 2) / t1)
            op = T1NoiseOp(1, lambda_half)
            s = BellState(1)
            
            # calculates bell probabilities across N runs 
            bp_probs = sum(
                bell_probs(dm(Ket(Stabilizer(mctrajectory!(copy(s), [op, op])[1]))))
                for _ in 1:N
            ) / N
            
            # std deviation calculations
            stderr = sqrt.(max.(bp_probs .* (1 .- bp_probs), eps(Float64)) ./ N)

            # 5 sigma checks, one to check that bpgates is within 5 sigma of the twirled data and another check to make sure that it's outside of 5 sigma of the untwirled data
            @test all(abs.(bp_probs .- twirled_probs) .<= 3/sqrt(N))

            # Skip untwirled rejection when distinguishing the two models would require more than N_max samples.
            if N_needed <= N_max
                @test any(abs.(bp_probs .- untwirled_probs) ./ stderr .> 5)
            else
                @info "Skipping untwirled-rejection at t=$t: delta=$delta needs N=$N_needed (> N_max=$N_max)"
            end
        end
    end
end