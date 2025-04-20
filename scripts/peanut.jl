using ProgressMeter
import FVADE
using Plots, LaTeXStrings


m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    V=x -> sum(x .^ 2) / 2,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> (1 - s)
)

h = 2^-3
is_in_Omega(x) = ((x[1]^2 - 3.9)^2 + x[2]^2 < 4^2)
limits = [(-20, 20), (-20, 20)]
println("Meshing")
mesh = FVADE.MeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)

println("size Ih = ", length(mesh.Ih))

# ρ0(x) = 1.0 * ((x[1] - 2.0)^2 + x[2]^2 < 2) + 1.0 * ((x[1] + 2.0)^2 + x[2]^2 < 2)
ρ0(x) = 0.6

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
length(ρ)
τ = h^2
@show τ

T = 2.0
N = ceil(Int64, T / τ)

println("Solving")
p = Progress(N)
M = maximum(ρ)
anim = @animate for n in 1:N
    FVADE.plot_2d(ρ, mesh)
    zlims!(0, 1.1)
    title!("t=$((n-1)*τ)")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )

    next!(p)
end
theplot = FVADE.plot_2d(ρ, mesh)
frame(anim, theplot)

mp4(anim, "figures/pme-peanut-saturation.mp4")

thetitle = L"U=s^2, V=|x|^2/2, K=0, \rho_0 = 0.6" * "\n" * latexstring("h=$h, τ=h^2, t=$T")
theplot = FVADE.plot_2d(ρ, mesh)
title!(thetitle)
savefig("figures/pme-peanut-saturation.pdf")