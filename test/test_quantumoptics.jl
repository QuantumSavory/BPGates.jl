@testitem "QuantumOptics comparisons" begin

using BPGates, QuantumClifford, QuantumOptics

using BPGates: T1NoiseOp

using LinearAlgebra: diag

# define some BP objects
λ = 0.2
N = 100000
s = BellState(1)
opBP = T1NoiseOp(1, λ)

# compute using BP the density matrix after T1 noise
ρBP = sum([dm(Ket(Stabilizer(apply!(s,op)))) for i in 1:N])/N

# define some QO objects
b = SpinBasis(1//2)
l0 = spinup(b)
l1 = spindown(b)
l00 = l0⊗l0
l01 = l0⊗l1
l10 = l1⊗l0
l11 = l1⊗l1
bell00 = (l00+l11)/sqrt(2)
bell01 = (l01+l10)/sqrt(2)
bell10 = (l00-l11)/sqrt(2)
bell11 = (l01-l10)/sqrt(2)
#|`00`|`+XX +ZZ`|`∣00⟩+∣11⟩`|`∣++⟩+∣--⟩`|`∣i₊i₋⟩+∣i₋i₊⟩`|
#|`01`|`+XX -ZZ`|`∣01⟩+∣10⟩`|`∣++⟩-∣--⟩`|`∣i₊i₊⟩-∣i₋i₋⟩`|
#|`10`|`-XX +ZZ`|`∣00⟩-∣11⟩`|`∣+-⟩+∣-+⟩`|`∣i₊i₊⟩+∣i₋i₋⟩`|
#|`11`|`-XX -ZZ`|`∣01⟩-∣10⟩`|`∣+-⟩-∣-+⟩`|`∣i₊i₋⟩-∣i₋i₊⟩`|

# compute using QO the density matrix after T1 noise
krausOp1 = projector(l0) + √(1-λ) * projector(l1)
krausOp2 = √(λ) * projector(l0, l1')
id = identityoperator(b)
k1 = krausOp1⊗id
k2 = krausOp2⊗id
k3 = id⊗krausOp1
k4 = id⊗krausOp2
ρ0 = dm(bell00)
ρ1  = k1*ρ0*k1' + k2*ρ0*k2'
ρQO = k3*ρ1*k3' + k4*ρ1*k4'

# switch to the Bell basis
to_bell_basis = projector(l00,bell00')+projector(l01,bell01')+projector(l10,bell10')+projector(l11,bell11')

ρbBP = to_bell_basis*ρBP*to_bell_basis'
ρbQO = to_bell_basis*ρQO*to_bell_basis'

@test_broken isapprox(diag(ρbBP.data), diag(ρbQO.data), atol=10/sqrt(N))

end
