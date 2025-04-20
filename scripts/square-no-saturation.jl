using ProgressMeter
import FVADE
using Plots

h = 2^-2
is_in_Omega(x) = (maximum(abs.(x)) < 4)
limits = [(-20, 20), (-20, 20)]
m = 2
U(s) = s^m / (m - 1)
Uprime(s) = m / (m - 1) * s^(m - 1)
V(x) = sum(x .^ 2) / 2
problem = FVADE.ADEProblem(
    U=U,
    Uprime=Uprime,
    V=V,
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

# ENV["JULIA_DEBUG"] = all

# ρ0(x) = 1.0 * (sum(x .^ 2) < 0.5)

# ρ0(x) = max((m - 1) / m * (1.0 - V(x)), 0.0)^(m - 1)

σ = 2^-1
ρ0(x) = 4.0 / σ^(FVADE.dimension(mesh)) * exp(-sum(x .^ 2) / σ)

# ρ0(x) = 1.0

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
τ = min(2^-6, h^2)
@show τ

T = 2.0
N = ceil(Int64, T / τ)

# N = 30

println("Solving")
p = Progress(N)
M = maximum(ρ)
anim = @animate for n in 1:N
    FVADE.plot_2d(ρ, mesh)
    zlims!(0, M)
    title!("t=$((n-1)*τ)")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
    # println("\n sum(ρ)=", sum(ρ))

    next!(p)
end

mp4(anim, "figures/pme-square.mp4")
