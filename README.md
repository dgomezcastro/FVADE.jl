# `FVADE.jl`

This packages implements a finite-volume scheme for aggregation-diffusion equations.

## Installation

Clone this repository to a folder, open julia in that folder and to install all dependencies run
```julia
using Pkg; 
Pkg.activate(".")
Pkg.instantiate()
```

## Simple examples of use
The `examples` folder contains simple examples of use $d=1$ and $d=2$. In you are running julia from the main folder to run the $d=1$ run 
```julia
include("examples/d=1.jl")
```

## Reproducing the figures in the paper
For an example run
```julia
include("gc2025/figure1-steady-state.jl")
include("gc2025/figure2-energy-decay.jl")
include("gc2025/figure3-convergence-Barenblatt.jl")
include("gc2025/figure4-convergence-halfh.jl")
```