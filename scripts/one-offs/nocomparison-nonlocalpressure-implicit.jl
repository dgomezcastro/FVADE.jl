using ProgressMeter
import FVADE
using Plots

using LinearAlgebra

h = 2^-2
@show h

is_in_Omega(x) = (abs(x[1]) < 7 && abs(x[2]) < 7)
limits = [(-10, 10), (-10, 10)]
m = 3
d = length(limits)
s = 0.5
problem = FVADE.ADEProblem(
    K=(x, y) -> (max(h, norm(x - y)))^(-2 + 2 * s),
    mobup=ρ -> ρ^(m - 1)
)
println("Meshing")
mesh = FVADE.UniformMeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)

println("size Ih = ", length(mesh.Ih))

# ENV["JULIA_DEBUG"] = all
bump(x, x0, r) = max(r^2 - norm(x - x0)^2, 0.0) / r^2

ρ01(x) = bump(x, [1.0, 0.0], 1.0)
ρ02(x) =
    ρ01(x) +
    5.0 * bump(x, [-2.0, 0.0], 1.0)

ρ1 = [Float64(ρ01(FVADE.x(i, h))) for i in mesh.Ih]
ρ2 = [Float64(ρ02(FVADE.x(i, h))) for i in mesh.Ih]
τ = h^2
@show τ

T = 0.5
N = ceil(Int64, T / τ)

println("Solving")
p = Progress(N)
M = max(maximum(ρ1), maximum(ρ2))
anim = @animate for n in 1:N+1
    p1 = FVADE.plot_2d(ρ1, mesh)
    zlims!(p1, 0, M)
    zlabel!(p1, "ρ1")
    p2 = FVADE.plot_2d(ρ2, mesh)
    zlims!(p2, 0, M)
    zlabel!(p2, "ρ2")
    p3 = FVADE.plot_2d((ρ2 - ρ1) / maximum(abs.(ρ2 - ρ1)), mesh, alpha=0.75)
    zlims!(p3, -M, M)
    zlabel!(p3, "(ρ2 - ρ1)/||ρ2-ρ1||_∞")

    p4 = plot(p1, p2, p3,
        plot_title="t=$(round((n-1)*τ,digits=3))",
        layout=(2, 2),
        size=(802, 608))

    global ρ1 = FVADE.iterate(
        ρ1, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
    global ρ2 = FVADE.iterate(
        ρ2, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )

    next!(p)
end

mp4(anim, "figures/nonlocalpressure-nocomparison.mp4")
