JET_flag = false

if get(ENV, "JET_TEST", "") == "true"
    JET_flag = true
else
    @info "Skipping JET tests -- must be explicitly enabled."
end

using Pkg
JET_flag && Pkg.add("JET")

using BPGates
using TestItemRunner

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]

  if JET_flag
    return :jet in ti.tags
  else
    push!(exclude, :jet)
  end

  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
