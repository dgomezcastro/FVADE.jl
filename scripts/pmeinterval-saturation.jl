using ProgressMeter
import FVADE
using Plots

h = 2^-4
is_in_Omega(x) = (maximum(abs.(x)) < 4)
limits = [(-20, 20)]
m = 2
U(s) = s^m / (m - 1)
Uprime(s) = m / (m - 1) * s^(m - 1)
V(x) = sum(x .^ 2)
problem = FVADE.ADEProblem(
    U=U,
    Uprime=Uprime,
    V=V,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1 - s
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

ρ0(x) = 0.6

ρ = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]
τ = h^2
@show τ

T = 5.0
N = ceil(Int64, T / τ)

# N = 30

println("Solving")
p = Progress(N)
M = 1.0
anim = @animate for n in 1:N+1
    plot(h * [i[1] for i in mesh.Ih], ρ, label="", linewidth=2)
    ylims!(0, M)
    title!("t=$(round((n-1)*τ,digits=3))")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
    # println("\n sum(ρ)=", sum(ρ))

    next!(p)
end

mp4(anim, "figures/pme-interval-saturation.mp4")
