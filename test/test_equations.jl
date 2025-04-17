@testset "ξ, v, F, L, iterate with trivial data" begin
    h = 1
    is_in_Omega(x) = (sum(x .^ 2) < 1.5)
    limits = [(-2, 2), (-2, 2)]
    problem = FVADE.ADEProblem(U=s -> s^2, V=x -> 0.0, K=(x, y) -> 0.0; Uprime=s -> 2 * s)
    mesh = FVADE.MeshADE(problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits)
    ρ = zeros(length(mesh.Ih))
    ξ = FVADE.ξ(ρ, problem, mesh)
    @test ξ == zeros(length(mesh.Ih))
    v = FVADE.v(ξ, mesh)
    @test v == zeros(length(mesh.Ih), FVADE.dimension(mesh))
    F = FVADE.F(ρ, v, problem, mesh)
    @test F == zeros(length(mesh.Ih), FVADE.dimension(mesh))
    τ = 1e-2
    L = FVADE.L(F, mesh, τ)
    @test L == zeros(length(mesh.Ih))
    ρ_next = FVADE.iterate(ρ, problem, mesh, τ)
    @test ρ_next == zeros(length(mesh.Ih))
end

@testset "ξ, v in d=1" begin
    h = 2^-2
    is_in_Omega(x) = minimum(-1 .< x .< 1)
    limits = [(-1.0, 1.0)]
    problem = FVADE.ADEProblem(U=s -> 0.0, V=x -> sum(x .^ 2) / 2, K=nothing, Uprime=s -> 0.0)
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
    F = FVADE.F(ρ, v, problem, mesh)
    @test F ≈ [0.0; 0.0; 0.0; -0.0; -0.0; -0.0; 0.0;;]
    τ = 0.1
    L = FVADE.L(F, mesh, τ)
    @test L ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end