push!(LOAD_PATH,"../src/")

using Documenter
#using DocumenterCitations
using BPGates

DocMeta.setdocmeta!(BPGates, :DocTestSetup, :(using QuantumClifford, BPGates); recursive=true)

#bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

makedocs(
#bib,
doctest = false,
clean = true,
sitename = "BPGates.jl",
format = Documenter.HTML(),
modules = [BPGates],
authors = "Shu Ge, Stefan Krastanov",
pages = [
"BPGates.jl" => "index.md",
"API" => "API.md",
]
)

deploydocs(
    repo = "github.com/Krastanov/BPGates.jl.git"
)
