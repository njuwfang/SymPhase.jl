import stim
from time import time

circuit = stim.Circuit.from_file(f"./stim_benchmark/random10.stim")
sampler = circuit.compile_sampler()
sampler.sample(shots=10000, bit_packed=True)

for (dep,n_cnot) in [(False,False), (False, True), (True, True)]:
    print(f"dep = {dep}, CNOT = {n_cnot}")
    if dep and n_cnot:
        suffix1 = "CNOT_dep"
    elif n_cnot:
        suffix1 = "CNOT"
    elif dep:
        suffix1 = "_dep"
    else:
        suffix1 = ""

    with open(f"./benchmark_stim{suffix1}.dat", "w") as io:
        print("nq init_time sample_time", file=io)
        for j in range(10,1001,10):
            if n_cnot:
                suffix2 = f"_{j>>1}"
            else:
                suffix2 = ""

            t0 = time()
            circuit = stim.Circuit.from_file(f"./stim_benchmark/random{j}{suffix2}{suffix1}.stim")
            sampler = circuit.compile_sampler()
            t1 = time()
            sampler.sample(shots=10000, bit_packed=True)
            t2 = time()
            print(f"{j} {t1-t0} {t2-t1}", file=io)
            print(f"{j} {t1-t0} {t2-t1}")

