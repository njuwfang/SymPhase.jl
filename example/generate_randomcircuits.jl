ENV["JULIA_PKG_USING_AUTOINSTALL"] = "yes"
using StatsBase: sample

mkpath("./stim_benchmark/")

function generate_random_circuit(n;depolarization=false,n_cnot=false)
    circuit = ""
    for j in 1:n
        circuit *= "H$(join(" $(k)" for k in sample(1:n,div(n,3)+1,replace=false)))\n"
        circuit *= "S$(join(" $(k)" for k in sample(1:n,div(n,3)+1,replace=false)))\n"
        if depolarization
          circuit *= "DEPOLARIZE1(0.01)$(join(" $(k)" for k in 1:n))\n"
        end
        if n_cnot
          circuit *= "CX$(join(" $(k)" for k in sample(1:n,n,replace=false)))\n"
        else
          circuit *= "CX$(join(" $(k)" for k in sample(1:n,10,replace=false)))\n"
        end
        circuit *= "M$(join(" $(k)" for k in sample(1:n,div(n,20)+1,replace=false)))\n"
    end
    circuit *= "M$(join(" $(k)" for k in 1:n))"
    circuit
end

for (dep,n_cnot) in [(false,false), (false, true), (true, true)]
@info "dep = $dep, CNOT = $n_cnot"
for j in 10:10:1000
    @info "j = $j"
    suffix = dep ? n_cnot ? "_$(j>>1)CNOT_dep" : "_dep" : n_cnot ? "_$(j>>1)CNOT" : ""
    a = generate_random_circuit(j;depolarization=dep,n_cnot=n_cnot)
    filename = "./stim_benchmark/random$(j)$(suffix).stim"
    open(filename, "w") do io
        println(io, a)
    end
end
end
