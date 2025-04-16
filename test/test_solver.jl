@testset "Newton solver" begin
    f(x) = x .^ 2
    y = FVADE.Newton(f, [0.0])
    @test y ≈ [0.0]

    y = FVADE.Newton(f, [1e-2]; abs_tol=1e-6)
    @test abs(y[1]) < 1e-5
end