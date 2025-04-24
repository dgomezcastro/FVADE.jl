import FVADE
using ProgressMeter
using Plots, LaTeXStrings, Plots.Measures
using LinearAlgebra
using JLD2

# ρ0(x) = 0.6
# ρ0_text = latexstring("$(ρ0(0.0))")

ρ0(x) = min(2.0 * exp(-norm(x / 0.5)^2), 1.0)
ρ0_text = L"min\{ 2 e^{-|x/0.5|^2}, 1 \}"

## d=1
# is_in_Omega(x) = (-4 < x[1] < 4)
# limits = [(-20, 20)]

## d=2
is_in_Omega(x) = (-1 < x[1] < 1 && -1 < x[2] < 1)
limits = [(-20, 20), (-20, 20)]

exponent_of_tau::Integer = 2
hs = 2.0 .^ collect(-2:-1:-5)

d = length(limits)

m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    # V=x -> sum(x .^ 2),
    V=nothing,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> (1 - s)
)
plottitle = latexstring("\\mathrm{m} = \\rho(1-\\rho), U=\\rho^2, V = 0, K = 0, \\rho_0 = ") * ρ0_text *
            "\n" * latexstring("\\Omega = [-4,4]^$d, T = 1, τ = h^{$exponent_of_tau}")

T = 0.5

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

L1errors = Vector(undef, length(hs))
Ihhalf, ρhhalf = solve(hs[1])
for (k, h) in enumerate(hs)
    @show ((k - 1) / (length(hs) - 1))
    Ih, ρh = Ihhalf, ρhhalf
    τ = h^exponent_of_tau

    global Ihhalf, ρhhalf = solve(h / 2)

    meshmapper = [findfirst(==(2 * i), Ihhalf) for i in Ih]

    L1errors[k] = τ * h^d * sum(abs.(ρhhalf[1:2^exponent_of_tau:end, meshmapper] - ρh))
    @show h, L1errors[k]
end


plot(hs, L1errors,
    scale=:log10,
    marker=:circle,
    xlabel=L"h",
    ylabel=L"ε_h^{(2)}",
    label="Simulation",
    title=plottitle,
    size=(700, 400),
    topmargin=10mm,
    bottommargin=5mm,
    leftmargin=5mm
)
plot!(hs, hs, label=L"Linear scaling $\varepsilon_h = h$")

savefig("figures/convergence-d$d.pdf")

save("figures/convergence-d$d.jld2", Dict("title" => plottitle, "hs" => hs, "L1errors" => L1errors))