import FVADE
using ProgressMeter
using Plots, LaTeXStrings, Plots.Measures
using LinearAlgebra
using JLD2

L = 4

d = 2

if d == 1
    is_in_Omega(x) = (-L < x[1] < L)
    limits = [(-L, L)]
    hs = 2.0 .^ (-2:-0.5:-5)
elseif d == 2
    is_in_Omega(x) = (-L < x[1] < L && -L < x[2] < L)
    limits = [(-L, L), (-L, L)]
    hs = 2.0 .^ (-2:-0.5:-4)
end

title = "convergence-Barenblatt-$(d)d"
m = 2
mass = 2.0
B = FVADE.BarenblattPME(m=Float64(m), dimension=d, mass=mass)

t0 = 1.0
ρ0(x) = B(x, t0)
ρ_exact(x, t) = B(x, t + t0)
ρ0_text = latexstring("\$\\textrm{Barenblatt}(M=$mass, t_0=$t0)\$")

exponent_of_tau::Integer = 2

T = 1.0

problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    # V=x -> sum(x .^ 2),
    V=nothing,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1
)



function solve(h)
    mesh = FVADE.UniformMeshADE(
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    τ = h^exponent_of_tau
    t = 0:τ:T
    ρ = Matrix{Float64}(undef, length(t), length(mesh.Ih))

    ρ[1, :] = FVADE.initialize_ρ(ρ0, mesh)


    @showprogress for n in 2:length(t)
        global ρ[n, :] = FVADE.iterate(
            ρ[n-1, :],
            problem, mesh, τ; abs_tol=1e-8, max_iters=10
        )
    end
    return mesh.Ih, t, ρ
end

L1errors = Vector(undef, length(hs))
for (k, h) in enumerate(hs)
    @show ((k - 1) / (length(hs) - 1))
    τ = h^exponent_of_tau

    global Ih, t, ρh = solve(h)

    ρ_exact_on_Ih = [ρ_exact(FVADE.x(i, h), t) for t in t, i in Ih]

    τ = t[2] - t[1]
    L1errors[k] = τ * h^d * sum(abs.(ρ_exact_on_Ih - ρh))
    @show h, L1errors[k]
end

plottitle = latexstring("\\mathrm{m} = \\rho, U=\\rho^$m, V = 0, K = 0, \\rho_0 = ") * ρ0_text *
            "\n" * latexstring("\\Omega = [-$L,$L]^$(d), T = $T, τ = h^{$exponent_of_tau}")

plot(hs, L1errors,
    scale=:log10,
    marker=:circle,
    xlabel=L"h",
    ylabel=L"\varepsilon^{(1)}_h",
    label="Simulation",
    title=plottitle,
    size=(700, 400),
    topmargin=10mm,
    bottommargin=5mm,
    leftmargin=5mm
)
plot!(hs, hs, label=L"Linear scaling $\varepsilon_h = h$")

savefig("figures/$title.pdf")

save("figures/$title.jld2", Dict("plottitle" => plottitle, "hs" => hs, "L1errors" => L1errors))