@testitem "QuantumOptics comparisons" begin

using Test

using BPGates, QuantumClifford, QuantumOpticsBase

using BPGates: T1NoiseOp, T2NoiseOp

using LinearAlgebra: diag

@testset "basis states" begin

# define some QO objects
b = SpinBasis(1//2)
l0 = spinup(b)
l1 = spindown(b)
l00 = l0⊗l0
l01 = l1⊗l0 # XXX be careful with the reversed endiandness - this is ket [0 1 0 0]
l10 = l0⊗l1 # XXX be careful with the reversed endiandness - this is ket [0 0 1 0]
l11 = l1⊗l1
bell00 = (l00+l11)/sqrt(2)
bell10 = (l00-l11)/sqrt(2) # XXX                                # bellstateindex = 2
bell01 = (l01+l10)/sqrt(2) # XXX                                # bellstateindex = 3
bell11 = (l01-l10)/sqrt(2)
#|`00`|`+XX +ZZ`|`∣00⟩+∣11⟩`|`∣++⟩+∣--⟩`|`∣i₊i₋⟩+∣i₋i₊⟩`|
#|`10`|`-XX +ZZ`|`∣00⟩-∣11⟩`|`∣+-⟩+∣-+⟩`|`∣i₊i₊⟩+∣i₋i₋⟩`|       # be careful : bellstateindex = 2
#|`01`|`+XX -ZZ`|`∣01⟩+∣10⟩`|`∣++⟩-∣--⟩`|`∣i₊i₊⟩-∣i₋i₋⟩`|       # be careful : bellstateindex = 3
#|`11`|`-XX -ZZ`|`∣01⟩-∣10⟩`|`∣+-⟩-∣-+⟩`|`∣i₊i₋⟩-∣i₋i₊⟩`|

@test BellState([0,0]) == BellState(BPGates.int_to_bit(1,2))
@test Ket(Stabilizer(BellState([0,0]))) ≈ Ket(S"XX ZZ") ≈ bell00

@test BellState([1,0]) == BellState(BPGates.int_to_bit(2,2))
@test Ket(Stabilizer(BellState([1,0]))) ≈ Ket(S"-XX ZZ") ≈ bell10

@test BellState([0,1]) == BellState(BPGates.int_to_bit(3,2))
@test Ket(Stabilizer(BellState([0,1]))) ≈ Ket(S"XX -ZZ") ≈ bell01

@test BellState([1,1]) == BellState(BPGates.int_to_bit(4,2))
@test Ket(Stabilizer(BellState([1,1]))) ≈ Ket(S"-XX -ZZ") ≈ -bell11

end

@testset "T1 and T2 noise" begin

# define some BP objects
λ = 0.2
N = 100000
opBP_T1 = T1NoiseOp(1, λ)
opBP_T2 = T2NoiseOp(1, λ)

# define some QO objects
b = SpinBasis(1//2)
l0 = spinup(b)
l1 = spindown(b)
l00 = l0⊗l0
l01 = l1⊗l0 # XXX be careful with the reversed endiandness - this is ket [0 1 0 0]
l10 = l0⊗l1 # XXX be careful with the reversed endiandness - this is ket [0 0 1 0]
l11 = l1⊗l1
bell00 = (l00+l11)/sqrt(2)
bell10 = (l00-l11)/sqrt(2) # XXX                                # bellstateindex = 2
bell01 = (l01+l10)/sqrt(2) # XXX                                # bellstateindex = 3
bell11 = (l01-l10)/sqrt(2)
#|`00`|`+XX +ZZ`|`∣00⟩+∣11⟩`|`∣++⟩+∣--⟩`|`∣i₊i₋⟩+∣i₋i₊⟩`|
#|`10`|`-XX +ZZ`|`∣00⟩-∣11⟩`|`∣+-⟩+∣-+⟩`|`∣i₊i₊⟩+∣i₋i₋⟩`|       # be careful : bellstateindex = 2
#|`01`|`+XX -ZZ`|`∣01⟩+∣10⟩`|`∣++⟩-∣--⟩`|`∣i₊i₊⟩-∣i₋i₋⟩`|       # be careful : bellstateindex = 3
#|`11`|`-XX -ZZ`|`∣01⟩-∣10⟩`|`∣+-⟩-∣-+⟩`|`∣i₊i₋⟩-∣i₋i₊⟩`|

# T1 noise in the QO formalism
krausOp1_T1 = projector(l0) + √(1-λ) * projector(l1)
krausOp2_T1 = √(λ) * projector(l0, l1')
id = identityoperator(b)
k1_T1 = krausOp1_T1⊗id
k2_T1 = krausOp2_T1⊗id
k3_T1 = id⊗krausOp1_T1
k4_T1 = id⊗krausOp2_T1

# T2 noise in the QO formalism
krausOp1_T2 = √(1-λ/2) * (projector(l0)+projector(l1))
krausOp2_T2 = √(λ/2) * (projector(l0)-projector(l1))
id = identityoperator(b)
k1_T2 = krausOp1_T2⊗id
k2_T2 = krausOp2_T2⊗id
k3_T2 = id⊗krausOp1_T2
k4_T2 = id⊗krausOp2_T2

# switch to the Bell basis
to_bell_basis = projector(l00,bell00')+projector(l01,bell01')+projector(l10,bell10')+projector(l11,bell11')

for bellstateindex in 1:4

    # compute using BP the density matrix after T1 noise
    s = BellState(BPGates.int_to_bit(bellstateindex, Val(2)))
    ρBP_T1 = sum([dm(Ket(Stabilizer(apply!(copy(s),opBP_T1)))) for i in 1:N])/N
    ρBP_T2 = sum([dm(Ket(Stabilizer(apply!(copy(s),opBP_T2)))) for i in 1:N])/N

    # compute using QO the density matrix after T1 noise
    ψ = [bell00,bell10,bell01,bell11][bellstateindex]

    @test abs(ψ' * Ket(Stabilizer(s))) ≈ 1
    ρ0 = dm(ψ)
    ρ1_T1  = k1_T1*ρ0*k1_T1' + k2_T1*ρ0*k2_T1'
    ρQO_T1 = k3_T1*ρ1_T1*k3_T1' + k4_T1*ρ1_T1*k4_T1'

    ρbBP_T1 = to_bell_basis*ρBP_T1*to_bell_basis'
    ρbQO_T1 = to_bell_basis*ρQO_T1*to_bell_basis'

    ρ1_T2  = k1_T2*ρ0*k1_T2' + k2_T2*ρ0*k2_T2'
    ρQO_T2 = k3_T2*ρ1_T2*k3_T2' + k4_T2*ρ1_T2*k4_T2'

    ρbBP_T2 = to_bell_basis*ρBP_T2*to_bell_basis'
    ρbQO_T2 = to_bell_basis*ρQO_T2*to_bell_basis'

    @test isapprox(diag(ρbBP_T1.data), diag(ρbQO_T1.data), atol=10/sqrt(N))
    @test isapprox(diag(ρbBP_T2.data), diag(ρbQO_T2.data), atol=10/sqrt(N))

end
end

end
