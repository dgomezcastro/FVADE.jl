@testset "PME on circle" begin
    h = 0.5
    is_in_Omega(x) = (sum(x .^ 2) < 1.5)
    limits = [(-4, 4), (-4, 4)]
    V(x) = sum(x .^ 2) / 2
    problem = FVADE.ADEProblem(
        s -> s^2 / 2,
        V,
        (x, y) -> 0.0
    )
    mesh = FVADE.MeshADE(
        problem=problem,
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    ρ = ones(length(mesh.Ih))
    τ = 1e-1
    for _ in 1:5/τ
        ρ = FVADE.iterate(
            ρ, problem, mesh, τ
        )
    end

    ρ_exact(C, x) = max(C - V(x), 0.0)
    xs = [FVADE.x(i, mesh.h) for i in mesh.Ih]
    min_error = Inf
    C_min = 5
    for C in C_min:-0.1:0.0
        error = maximum(abs.(ρ - ρ_exact.(C, xs)))
        if error < min_error
            C_min = C
            min_error = error
        end
    end
    @test min_error < 5e-2
end