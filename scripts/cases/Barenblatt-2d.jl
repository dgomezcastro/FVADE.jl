L = 4

d = 2
is_in_Omega(x) = (-L < x[1] < L && -L < x[2] < L)
limits = [(-L, L), (-L, L)]

m = 2
mass = 2.0
B = FVADE.BarenblattPME(m=Float64(m), dimension=d, mass=mass)

t0 = 1.0
ρ0(x) = B(x, t0)
ρ_exact(x, t) = B(x, t + t0)
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