using Test, Random
using QuantumClifford
using BPGates

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

doset("quantumclifford")    && include("./test_quantumclifford.jl")
doset("jet")                && haskey(ENV,"QUANTUMCLIFFORD_JET_TEST") && ENV["QUANTUMCLIFFORD_JET_TEST"]=="true" && include("./test_jet.jl")
#TODO doset("allocations")        && VERSION >= v"1.7" && include("./test_allocations.jl")
doset("doctests")           && VERSION == v"1.7" && include("./doctests.jl")

using Aqua
doset("aqua") && begin
    Aqua.test_all(BPGates, ambiguities=false)
    Aqua.test_ambiguities([BPGates,Core]) # otherwise Base causes false positives
end
