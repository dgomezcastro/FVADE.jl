using ProgressMeter
import FVADE
using Plots

ρ0(x) = 0.6
is_in_Omega(x) = (-4 < x[1] < 4 && -4 < x[2] < 4)

h = 2^-3
limits = [(-20, 20), (-20, 20)]
m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    V=x -> sum(x .^ 2),
    K=nothing,
    mobup=s -> s,
    mobdown=s -> (1 - s)
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


# ρ0(x) = max((m - 1) / m * (1.0 - V(x)), 0.0)^(m - 1)

# σ = 2^-1
# ρ0(x) = min(4.0 / σ^(FVADE.dimension(mesh)) * exp(-sum(x .^ 2) / σ), 1.0)

# ρ0(x) = 1.0

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
τ = h^2
@show τ

T = 1.0
N = ceil(Int64, T / τ) + 1

# N = 30

println("Solving")
p = Progress(N)
M = 1.0
anim = @animate for n in 1:N
    FVADE.plot_2d(ρ, mesh)
    zlims!(0, M)
    title!("t=$(round((n-1)*τ,digits=3))")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
    # println("\n sum(ρ)=", sum(ρ))

    next!(p)
end
savefig("figures/pme-squre-satuaration.pdf")
mp4(anim, "figures/pme-square-saturation.mp4")
