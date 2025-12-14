push!(LOAD_PATH,"../src/")

using Documenter
using AnythingLLMDocs
#using DocumenterCitations
using BPGates

DocMeta.setdocmeta!(BPGates, :DocTestSetup, :(using QuantumClifford, BPGates); recursive=true)

doc_modules = [BPGates]

api_base="https://anythingllm.krastanov.org/api/v1"
anythingllm_assets = integrate_anythingllm(
    "BPGates",
    doc_modules,
    @__DIR__,
    api_base;
    repo = "github.com/QuantumSavory/BPGates.jl.git",
    options = EmbedOptions(),
)

#bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

makedocs(
#bib,
doctest = false,
clean = true,
warnonly = :missing_docs,
sitename = "BPGates.jl",
format = Documenter.HTML(assets = anythingllm_assets),
modules = doc_modules,
authors = "Shu Ge, Stefan Krastanov",
pages = [
"BPGates.jl" => "index.md",
"Most Common Gates" => "most_common.md",
"Coincidence Measurements" => "measurements.md",
"Noisy Gates" => "noisy.md",
"API" => "API.md",
]
)

deploydocs(
    repo = "github.com/QuantumSavory/BPGates.jl.git"
)
