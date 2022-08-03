using Documenter

function doctests()
    @testset "Doctests" begin
        DocMeta.setdocmeta!(BPGates, :DocTestSetup, :(using QuantumClifford, BPGates); recursive=true)
        doctest(BPGates)
    end
end

doctests()
