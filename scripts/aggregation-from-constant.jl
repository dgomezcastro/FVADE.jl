import FVADE
using ProgressMeter
using Plots, LaTeXStrings
using LinearAlgebra

W(x) = norm(x)^2 / 2 #- norm(x)^4 / 4

ρ0(x) = 0.6
is_in_Omega(x) = (-4 < x[1] < 4 && -4 < x[2] < 4)
h = 2^-2
exponent_of_tau = 3
τ = h^exponent_of_tau

limits = [(-20, 20), (-20, 20)]
d = length(limits)
m = 2
U = s -> s^m / (m - 1)
Uprime = s -> m / (m - 1) * s^(m - 1)
K = (x, y) -> W(x - y)
problem = FVADE.ADEProblem(
    U=U, Uprime=Uprime,
    V=nothing,
    K=K,
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

@show τ

T = 0.5
N = ceil(Int64, T / τ) + 1

# N = 30

println("Solving")
p = Progress(N)
M = 1.0

free_energies = Vector{Float64}(undef, N)
anim = @animate for n in 1:N
    free_energies[n] = FVADE.free_energy(ρ, problem, mesh)

    FVADE.plot_2d(ρ, mesh)
    zlims!(0, M)
    title!("t=$(round((n-1)*τ,digits=3))")

    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
    )
    # println("\n sum(ρ)=", sum(ρ))

    next!(p)
end
mp4(anim, "figures/pme-square-aggregation.mp4")

thetitle = latexstring("U=s^{$m}, V=0, K=|x-y|^2/2, \\rho_0 = 0.6") * "\n" * latexstring("h=$h, τ=h^$exponent_of_tau, t=$T")
title!(thetitle)
savefig("figures/pme-square-aggregation.pdf")

thetitle = latexstring("U=s^{$m}, V=0, K=|x-y|^2/2, \\rho_0 = 0.6") * "\n" * latexstring("\\Omega=[-4,4]^2, h=$h, τ=h^$exponent_of_tau")
plot(τ * [0:(N-1);], free_energies,
    title=thetitle,
    label="",
    xlabel=L"t",
    ylabel=L"\mathcal{F}[\rho_t]",
    linewidth=2)
savefig("figures/pme-square-aggregation-energy.pdf")
