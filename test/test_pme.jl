@testset "PME on circle with positive data" begin
    h = 0.5
    is_in_Omega(x) = (sum(x .^ 2) < 1.5)
    limits = [(-4, 4), (-4, 4)]
    V(x) = sum(x .^ 2) / 2
    problem = FVADE.ADEProblem(
        U=s -> s^2 / 2,
        Uprime=s -> s,
        V=V,
        K=nothing
    )
    mesh = FVADE.UniformMeshADE(
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    ρ = ones(length(mesh.Ih))
    τ = 1e-1
    for _ in 1:5/τ
        ρ = FVADE.iterate(
            ρ, problem, mesh, τ, abs_tol=1e-5
        )
    end

    ρ_exact(C, x) = max(C - V(x), 0.0)
    xs = [FVADE.x(i, mesh.h) for i in mesh.Ih]
    min_error = Inf
    C_min = 20
    for C in C_min:-1e-3:0.0
        error = maximum(abs.(ρ - ρ_exact.(C, xs)))
        if error < min_error
            C_min = C
            min_error = error
        end
    end
    @test min_error < 1e-3
end