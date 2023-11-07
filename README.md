# SymPhase.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://njuwfang.github.io/SymPhase.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://njuwfang.github.io/SymPhase.jl/dev/)
[![Build Status](https://github.com/njuwfang/SymPhase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/njuwfang/SymPhase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/njuwfang/SymPhase.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/njuwfang/SymPhase.jl)

A high-performance Julia module for sampling stabilizer circuits by symbolizing the phases of Pauli strings.
See [the arXiv paper describing SymPhase](https://arxiv.org).

## Installation

To install it, run [Julia](https://julialang.org/) REPL and use:

```
] add https://github.com/njuwfang/SymPhase.jl.git
```

## Benchmarks

We have compared SymPhase.jl with Stim in [the arXiv paper](https://arxiv.org).

To reproduce the experiment, follow this procedure

1. Install Stim (das of 7/11/2023, the latest stable version is 1.12.0) with `pip`.
    ```
    pip install stim==1.12.0
    ```
2. Go to the `example` directory and generate benchmark files in `.stim` format.
    ```
    julia generate_randomcircuits.jl
    ```

    Each circuit is made up of n qubits with n layers. Each layer randomly applies an H, S and I gate to each qubit, then applies CNOT gates, then samples 5\% of the qubits to measure in the computational basis. At the end of the circuit, each qubit is measured in the computational basis.
    - for `./stim_benchmark/randomxxx.stim`, each layer randomly selects 5 pairs (1 pair in Stim's paper) of qubits to apply CNOT gates.
    - for `./stim_benchmark/randomxxx_xxxCNOT.stim`, each layer randomly selects $\lfloor\frac{n}{2}\rfloor$ pairs of qubits to apply CNOT gates.
    - for `./stim_benchmark/randomxxx_xxxCNOT_dep.stim`, each layer randomly selects $\lfloor\frac{n}{2}\rfloor$ pairs of qubits to apply CNOT gates and additionally applies single-qubit depolarize noise to each qubit.

3. Benchmark Stim and SymPhase.jl.
    ```
    python benchmark_stim.py
    ```
    ```
    julia benchmark_symphase.jl
    ```
    The benchmark results are stored in files with the `.dat` suffix.