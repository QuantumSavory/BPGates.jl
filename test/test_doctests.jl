@testitem "doctests" begin

using Documenter
using BPGates

DocMeta.setdocmeta!(BPGates, :DocTestSetup, :(using QuantumClifford, BPGates); recursive=true)
doctest(BPGates)

end
