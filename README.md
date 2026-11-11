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

## Choosing the nonlinear solver
Each implicit time step `FVADE.iterate` solves a nonlinear system with [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl). The keyword `solver` chooses the method:
- `solver=:SimpleDFSane` (default): a derivative-free quasi-Newton method, which is significantly faster.
- `solver=:NewtonRaphson`: Newton-Raphson's method, with the Jacobian computed by `ForwardDiff.jl`, and at most `max_iters` iterations.

For example
```julia
ρ = FVADE.iterate(ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=10, solver=:NewtonRaphson)
```

## Reproducing the figures in the paper
The figures in the paper are computed with `solver=:NewtonRaphson`.

For an example run
```julia
include("gc2025/figure1-steady-state.jl")
include("gc2025/figure2-energy-decay.jl")
include("gc2025/figure3-convergence-Barenblatt.jl")
include("gc2025/figure4-convergence-halfh.jl")
```
The convergence scripts save the errors in `figures/*.jld2` and then plot them. To redo only the plots from the saved errors, without computing again, run
```julia
include("gc2025/figure3-plot.jl")
include("gc2025/figure4-plot.jl")
```