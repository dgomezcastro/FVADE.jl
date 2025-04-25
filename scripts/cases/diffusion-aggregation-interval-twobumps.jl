W(x) = norm(x)^2 / 2 #- norm(x)^4 / 4

bump(x, x0, r) = max(r^2 - norm(x - x0)^2, 0.0) / r^2

ρ0(x) = bump(x, [1.0], 0.5) + bump(x, [-1.0], 0.5)

is_in_Omega(x) = (-4 < x[1] < 4)
h = 2^-3
exponent_of_tau = 2
τ = h^exponent_of_tau

limits = [(-20, 20)]
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

plottitle = latexstring("\\mathrm{m}=\rho(1-\rho), U=s^{$m}, V=0, K=|x-y|^2/2") * "\n" * latexstring("\\Omega=[-$L,$L]^{$d}, h=$h, τ=h^$exponent_of_tau")
