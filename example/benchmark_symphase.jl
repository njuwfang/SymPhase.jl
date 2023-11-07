using Pkg
Pkg.activate("../")

using SymPhase


c = parse_stim("./stim_benchmark/random10.stim")
sampler = Sampler(c)

for (dep,n_cnot) in [(false,false), (false, true), (true, true)]
  @info "dep = $dep, CNOT = $n_cnot"
  suffix1 = dep ? n_cnot ? "CNOT_dep" : "_dep" : n_cnot ? "CNOT" : ""
open("benchmark_symphase$(suffix1).dat", "w") do io
  println(io, "nq init_time sample_time")
  for j in 10:10:1000
    suffix2 = n_cnot ? "_$(j>>1)" : ""
    t0 = time()
    c = parse_stim("./stim_benchmark/random$(j)$(suffix2)$(suffix1).stim")
    sampler = Sampler(c)
    t1 = time()
    sampler(10000)
    t2 = time()
    println(io, "$(j) $(t1-t0) $(t2-t1)")
    println("$(j) $(t1-t0) $(t2-t1)")
  end
end
end
