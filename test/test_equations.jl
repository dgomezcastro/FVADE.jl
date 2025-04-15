@testset "ξ, v, F, H, iterate with trivial data" begin
    h = 1
    is_in_Omega(x) = (sum(x .^ 2) < 1.5)
    limits = [(-2, 2), (-2, 2)]
    problem = FVADE.ADEProblem(s -> s^2, x -> 0.0, (x, y) -> 0.0)
    mesh = FVADE.MeshADE(problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits)
    ρ = zeros(length(mesh.Ih))
    ξ = FVADE.ξ(ρ, problem, mesh)
    @test ξ == zeros(length(mesh.Ih))
    v = FVADE.v(ξ, mesh)
    @test v == zeros(length(mesh.Ih), FVADE.dimension(mesh))
    F = FVADE.F(ρ, v, mesh)
    @test F == zeros(length(mesh.Ih), FVADE.dimension(mesh))
    τ = 1e-2
    H = FVADE.H(ρ, F, mesh, τ)
    @test H == zeros(length(mesh.Ih))
    ρ_next = FVADE.iterate(ρ, problem, mesh, τ)
    @test ρ_next == zeros(length(mesh.Ih))
end
