using ProgressMeter
import FVADE
using Plots, LaTeXStrings, Plots.Measures

ρ0(x) = 0.6

is_in_Omega(x) = (-4 < x[1] < 4 && -4 < x[2] < 4)
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

d = length(limits)
exponent_of_tau::Integer = 2

problem = FVADE.ADEProblem(
    U=s -> s^2,
    Uprime=s -> 2 * s,
    V=x -> sum(x .^ 2) / 2,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> (1 - s)
)
T = 1.0

function solve(h)
    mesh = FVADE.MeshADE(
        problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    τ = h^exponent_of_tau
    N = round(Int64, T / τ)
    ρ = Matrix{Float64}(undef, N, length(mesh.Ih))

    ρ[1, :] = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]

    @showprogress for n in 2:N
        global ρ[n, :] = FVADE.iterate(
            ρ[n-1, :],
            problem, mesh, τ; abs_tol=1e-8, max_iters=10
        )
    end
    return mesh.Ih, ρ
end

hs = 2.0 .^ collect(-2:-1:-5)
L1errors = Vector(undef, length(hs))
k = 1
@show (k / length(hs))
Ihhalf, ρhhalf = solve(hs[1])
for (k, h) in enumerate(hs)
    @show (k / length(hs))
    Ih, ρh = Ihhalf, ρhhalf
    τ = h^exponent_of_tau

    Ihhalf, ρhhalf = solve(h / 2)

    meshmapper = [findfirst(==(2 * i), Ihhalf) for i in Ih]

    L1errors[k] = τ * h^d * sum(abs.(ρhhalf[1:2^exponent_of_tau:end, meshmapper] - ρh))
end

plot(hs, L1errors,
    scale=:log10,
    marker=:circle,
    xlabel=L"h",
    ylabel=L"||ρ_h - ρ_{h / 2}||_{L^2( (0,T) \times \Omega )}",
    label="",
    title=latexstring("\\mathrm{m} = \\rho(1-\\rho), U=\\rho^2, V = |x|^2/2, K = 0, \\rho_0 = $(ρ_0([0,0]))") *
          "\n" * latexstring("\\Omega = [-4,4]^2, T = 1, τ = h^{$exponent_of_tau}"),
    size=(700, 300),
    topmargin=10mm
)

savefig("figures/convergence.pdf")