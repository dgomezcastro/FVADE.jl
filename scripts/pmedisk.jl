using ProgressMeter
import FVADE

h = 0.1
is_in_Omega(x) = (sum(x .^ 2) < 1.5)
limits = [(-4, 4), (-4, 4)]
V(x) = sum(x .^ 2) / 2
problem = FVADE.ADEProblem(
    s -> s^2 / 2,
    V,
    (x, y) -> 0.0
)
println("Meshing")
mesh = FVADE.MeshADE(
    problem=problem,
    is_in_Omega=is_in_Omega,
    h=h,
    mesh_limits=limits
)

println("Solving")
ρ0(x) = 1.0
ρ = [ρ0(FVADE.x(i, h)) for i in mesh.Ih]
τ = 1e-1
T = 2.0
@showprogress for _ in 1:T/τ
    global ρ = FVADE.iterate(
        ρ, problem, mesh, τ
    )
end

FVADE.plot_2d(ρ, mesh)