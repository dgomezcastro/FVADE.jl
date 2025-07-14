rL = 1

is_in_Omega(x) = (-L < x[1] < L && -L < x[2] < L)
limits = [(-L, L), (-L, L)]

ρ0(x) = min(2.0 * exp(-norm(x / 0.5)^2), 1.0)
ρ0_text = L"min\{ 2 e^{-|x/0.5|^2}, 1 \}"

h = 2^-3
exponent_of_tau::Integer = 2
τ = h^exponent_of_tau

T = 2.0

m = 2
problem = FVADE.ADEProblem(
    U=s -> s^m / (m - 1),
    Uprime=s -> m / (m - 1) * s^(m - 1),
    V=nothing,
    K=nothing,
    mobup=s -> s,
    mobdown=s -> (1 - s)
)

println("Meshing")
mesh = FVADE.UniformMeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)
println("size Ih = ", length(mesh.Ih))

zlimit = 1.0