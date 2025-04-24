import FVADE
using ProgressMeter
using Plots, LaTeXStrings, Plots.Measures
using LinearAlgebra
using JLD2

L = 4

## d=1
is_in_Omega(x) = (-L < x[1] < L)
limits = [(-L, L)]

d = length(limits)
title = "interval-pme-$(d)d"
m = 2
mass = 2.0
B = FVADE.BarenblattPME(m=Float64(m), dimension=d, mass=mass)

t0 = 1.0
ρ0(x) = FVADE.evaluate(B, x, t0)
ρ_exact(x::Vector, t) = FVADE.evaluate(B, x, t + t0)
ρ_exact(x::Number, t) = FVADE.evaluate(B, [x], t + t0)
ρ0_text = L"Barenblatt"

h = 2.0^-5
exponent_of_tau::Integer = 2

T = 1.0

m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    # V=x -> sum(x .^ 2),
    V=nothing,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1
)

plottitle = latexstring("\\mathrm{m} = \\rho, U=\\rho^$m, V = 0, K = 0, \\rho_0 = ") * ρ0_text *
            "\n" * latexstring("\\Omega = [-$L,$L]^$d, T = $T, τ = h^{$exponent_of_tau}")


function solve(h)
    mesh = FVADE.MeshADE(
        problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    τ = h^exponent_of_tau
    t = 0:τ:T
    ρ = Matrix{Float64}(undef, length(t), length(mesh.Ih))

    ρ[1, :] = [Float64(ρ0(FVADE.x(i, h))) for i in mesh.Ih]

    @showprogress for n in 2:length(t)
        global ρ[n, :] = FVADE.iterate(
            ρ[n-1, :],
            problem, mesh, τ; abs_tol=1e-8, max_iters=10
        )
    end
    return mesh.Ih, t, ρ
end

Ih, tt, ρh = solve(h)


xx = [FVADE.x(i, h)[1] for i in Ih]
anim = @animate for k = 1:10:length(tt)
    t = tt[k]
    plot(xx, ρh[k, :], label="Numerical solution")
    plot!(xx, ρ_exact.(xx, t), label="Exact solution", marker=:circle)
    # @show h * sum(ρh[k, :])
    # @show h * sum(ρ_exact.(xx, t))
    title!("t=$(round(t, digits=3))")
    ylims!(0, 1)
end

mp4(anim, "figures/Barenblatt-1d.mp4")