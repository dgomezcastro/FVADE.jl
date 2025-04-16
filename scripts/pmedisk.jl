using ProgressMeter
import FVADE
using Plots

h = 2^-2
is_in_Omega(x) = (sum(x .^ 2) < 1.5)
limits = [(-4, 4), (-4, 4)]
problem = FVADE.ADEProblem(
    U=s -> s^2 / 2,
    V=x -> sum(x .^ 2) / 2,
    K=nothing
)
println("Meshing")
mesh = FVADE.MeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)
println("size Ih = ", length(mesh.Ih))

ENV["JULIA_DEBUG"] = all

println("Solving")
# ρ0(x) = 0.25 * (sum(x .^ 2) < 0.5)
ρ0(x) = exp(-sum(x .^ 2) / 0.1)
# ρ0(x) = 1.0

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
τ = h^2
T = 2.0
N = ceil(Int64, T / τ)
p = Progress(N)
M = maximum(ρ)
anim = @animate for n in 1:N
    FVADE.plot_2d(ρ, mesh)
    zlims!(0, M)
    title!("t=$(round((n-1)*τ))")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=10
    )
    # println("\n sum(ρ)=", sum(ρ))
    next!(p)
end

mp4(anim, "figures/pme-disk.mp4")
