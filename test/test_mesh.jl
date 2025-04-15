
@testset "Generate Ih" begin
    @testset "Square" begin
        h = 2
        is_in_Omega(x) = (maximum(abs.(x)) < 1.5)
        limits = [(-1, 1), (-1, 1)]
        Ih = FVADE.generate_Ih(is_in_Omega, h, limits)
        @test Ih == [[0, 0]]
    end
    @testset "Circle" begin
        h = 1
        is_in_Omega(x) = (sum(x .^ 2) < 1.5)
        limits = [(-2, 2), (-2, 2)]
        Ih = FVADE.generate_Ih(is_in_Omega, h, limits)
        @test Ih == [[0, -1], [-1, 0], [0, 0], [1, 0], [0, 1]]
    end
end

@testset "Generate Mesh" begin
    @testset "Square" begin
        h = 2
        is_in_Omega(x) = (maximum(abs.(x)) < 1.5)
        limits = [(-1, 1), (-1, 1)]
        problem = FVADE.ADEProblem(s -> s^2, x -> 0.0, (x, y) -> 0.0)
        mesh = FVADE.MeshADE(
            problem=problem,
            is_in_Omega=is_in_Omega,
            h=h,
            mesh_limits=limits)
        @test mesh.h == h
        @test mesh.Ih == [[0, 0]]
        @test mesh.VV == [0.0]
        @test mesh.KK == [0.0;;]
    end
    @testset "Circle" begin
        h = 1
        is_in_Omega(x) = (sum(x .^ 2) < 1.5)
        limits = [(-2, 2), (-2, 2)]
        problem = FVADE.ADEProblem(s -> s^2, x -> 0.0, (x, y) -> 0.0)
        mesh = FVADE.MeshADE(
            problem=problem,
            is_in_Omega=is_in_Omega,
            h=h,
            mesh_limits=limits)
        @test mesh.h == h
        @test mesh.Ih == [[0, -1], [-1, 0], [0, 0], [1, 0], [0, 1]]
        @test mesh.neighbours_plus == [nothing 3; 3 nothing; 4 5; nothing nothing; nothing nothing]
        @test mesh.neighbours_minus == [nothing nothing; nothing nothing; 2 1; 3 nothing; nothing 3]
    end
end
