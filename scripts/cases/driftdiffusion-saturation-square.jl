h = 2^-3
L = 4
is_in_Omega(x) = (-L < x[1] < L && -L < x[2] < L)
limits = [(-L, L), (-L, L)]
d = length(limits)

exponent_of_tau = 2
τ = h^exponent_of_tau

m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    V=x -> sum(x .^ 2) / 2,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> 1 - s
)

zlimit = 1

ρ0(x) = 0.6

T = 2.0

plottitle = latexstring("\\mathrm{m}=\\rho(1-\\rho), U=s^{$m}, V=|x|^2/2, K=0, \\rho_0 = 0.6") * "\n" * latexstring("\\Omega=[-L,L]^{$d}, h=$h, τ=h^$exponent_of_tau")