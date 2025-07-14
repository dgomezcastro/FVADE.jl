using QuadGK
@testset "Barenblatt does not give an error in d=1,2" begin
    B = FVADE.BarenblattPME(m=2.0, dimension=1, mass=1.0)
    B(1.0, 1.0)
    B = FVADE.BarenblattPME(m=2.0, dimension=2, mass=1.0)
    B([1.0, 1.0], 1.0)
end
@testset "BarenblattPME has the correct mass" begin
    @testset "d=1 at t=1.0" begin
        for mass = 0.1:0.1:2.0
            B = FVADE.BarenblattPME(m=2.0, dimension=1, mass=mass)
            mass_approx = quadgk(x -> B(x, 1.0), -Inf, Inf, rtol=1e-10)[1]
            @test abs(mass_approx - mass) / mass < 1e-4
        end
    end
    @testset "d=1 at t=2.0" begin
        for mass = 0.1:0.1:2.0
            B = FVADE.BarenblattPME(m=2.0, dimension=1, mass=mass)
            mass_approx = quadgk(x -> B(x, 2.0), -Inf, Inf, rtol=1e-10)[1]
            @test abs(mass_approx - mass) / mass < 1e-4
        end
    end
end;