using SafeTestsets
using BPGates

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

macro doset(descr)
    quote
        if doset($descr)
            @safetestset $descr begin include("test_"*$descr*".jl") end
        end
    end
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@doset "quantumclifford"
@doset "quantikz"
get(ENV,"JET_TEST","")=="true" && @doset "jet"
@doset "doctests"
@doset "bpgates"

using Aqua
doset("aqua") && begin
    Aqua.test_all(BPGates)
end
