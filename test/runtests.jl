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
            @safetestset $descr begin 
            if $descr == "mixedstateop" || $descr == "paulizop"
                include("test_kraus_mem_noise_helpers.jl")
            end
            include("test_"*$descr*".jl") end
        end
    end
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

# include("test_kraus_mem_noise_helpers.jl")

@doset "quantumclifford"
@doset "quantikz"
get(ENV,"JET_TEST","")=="true" && @doset "jet"
@doset "doctests"
@doset "bpgates"
@doset "mixedstateop"
@doset "paulizop"

using Aqua
doset("aqua") && begin
    Aqua.test_all(BPGates)
end
