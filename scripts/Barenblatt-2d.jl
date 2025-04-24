using ProgressMeter
import FVADE
using Plots, LaTeXStrings
using LinearAlgebra

title = "Barenblatt-2d"

m = 2
B = FVADE.BarenblattPME(m=2.0, dimension=2, mass=3.0)

t0 = 1e-1
ρ0(x) = FVADE.evaluate(B, x, t0)
ρ_exact(x, t) = FVADE.evaluate(B, x, t + t0)
ρ0_text = L"Barenblatt"

is_in_Omega(x) = (-2 < x[1] < 2 && -2 < x[2] < 2)
limits = [(-20, 20), (-20, 20)]

h = 2^-4
exponent_of_tau::Integer = 1
τ = h^exponent_of_tau
@show h, τ

T = 1.0
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    V=nothing,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1
)
println("Meshing")
mesh = FVADE.MeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)
println("size Ih = ", length(mesh.Ih))

# ENV["JULIA_DEBUG"] = all

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
N = ceil(Int64, T / τ) + 1

# N = 30

println("Solving")
p = Progress(N)
M = 1.0
xx = -2:h:2
yy = -2:h:2

anim = @animate for n in 1:N
    t = (n - 1) * τ
    # surface(xx, yy, (x, y) -> ρ_exact([x, y], t))

    FVADE.plot_2d(ρ, mesh)
    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )

    zlims!(0, M)
    title!("t=$(round(t,digits=3))")
    next!(p)
end
mp4(anim, "figures/$title.mp4")

thetitle = L"U=s^2, V=0, K=0" * "\n" * latexstring("\rho_0=") * ρ0_text * latexstring(", h=$h, τ=h^$exponent_of_tau, t=$T")
theplot = FVADE.plot_2d(ρ, mesh)
title!(thetitle)
zlims!(0, M)
savefig("figures/$title.pdf")