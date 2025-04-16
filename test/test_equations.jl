@testset "ξ, v, F, H, iterate with trivial data" begin
    h = 1
    is_in_Omega(x) = (sum(x .^ 2) < 1.5)
    limits = [(-2, 2), (-2, 2)]
    problem = FVADE.ADEProblem(s -> s^2, s -> 2 * s, x -> 0.0, (x, y) -> 0.0)
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

@testset "ξ, v in d=1" begin
    h = 2^-2
    is_in_Omega(x) = minimum(-1 .< x .< 1)
    limits = [(-1.0, 1.0)]
    problem = FVADE.ADEProblem(s -> 0.0, s -> 0.0, x -> sum(x .^ 2) / 2, nothing)
    mesh = FVADE.MeshADE(problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits)
    @test mesh.Ih == [[-3], [-2], [-1], [0], [1], [2], [3]]
    ρ = zeros(length(mesh.Ih))
    ξ = FVADE.ξ(ρ, problem, mesh)
    @test ξ ≈ [0.28125, 0.125, 0.03125, 0.0, 0.03125, 0.125, 0.28125]
    v = FVADE.v(ξ, mesh)
    @test v ≈ [0.625; 0.375; 0.125; -0.125; -0.375; -0.625; 0.0;;]
    F = FVADE.F(ρ, v, mesh)
    @test F ≈ [0.0; 0.0; 0.0; -0.0; -0.0; -0.0; 0.0;;]
    τ = 0.1
    H = FVADE.H(ρ, F, mesh, τ)
    @test H ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end