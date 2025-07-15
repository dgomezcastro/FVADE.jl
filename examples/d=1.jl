using Pkg;
Pkg.activate(".");
using Plots
import FVADE

is_in_Omega(x) = (-4 < x[1] < 4)
h = 2^-2
τ = h^2
T = 5.0

ρ0(x) = 0.3

problem = FVADE.ADEProblem(
    U=s -> s^2,
    V=x -> sum(x .^ 2) / 2,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1.0,
)
mesh = FVADE.UniformMeshADE(
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=[(-4, 4)]
)

ρ = FVADE.initialize_ρ(ρ0, mesh)

@gif for n in 0:(T/τ)
    global ρ
    FVADE.plot_1d(ρ, mesh)
    title!("t=$(round(n*τ,digits=3))")
    ylims!(0, 1.1)

    ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
end