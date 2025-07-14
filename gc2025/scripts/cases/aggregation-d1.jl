L = 4
d = 1

is_in_Omega(x) = (-L < x[1] < L)
h = 2^-2
limits = [(-L, L)]

exponent_of_tau = 3
τ = h^exponent_of_tau

name = "energy-decay-aggregation-d$d"

ρ0(x) = 0.6

m = 2
U = s -> s^m / (m - 1)
Uprime = s -> m / (m - 1) * s^(m - 1)
W(x) = norm(x)^2 / 2 #- norm(x)^4 / 4
K = (x, y) -> W(x - y)
problem = FVADE.ADEProblem(
    U=U, Uprime=Uprime,
    V=nothing,
    K=K,
    mobup=s -> s,
    mobdown=s -> (1 - s)
)

T = 5.0

plottitle = latexstring("\\mathrm{m}=\\rho(1-\\rho), U=s^{$m}, V=0, K=|x-y|^2/2, \\rho_0 = 0.6") * "\n" * latexstring("\\Omega=[-$L,$L]^{$d}, h=$h, τ=h^$exponent_of_tau")
