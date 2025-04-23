using QuadGK, SpecialFunctions

# The mass of the Barenblatt can be computed explicitly (as in [Vazquez, 2007]) with the change of variable $y = C^{\frac 1 2} z$
# $$
# \begin{aligned}
#     \int_{\mathbb R} B_C(x,C) dx &=  \int_{\mathbb R} ( \tfrac{m}{m-1}( C - \tfrac{y^2}{2}  )  )_+^{\frac{1}{m-1}} dy \\
#     &= \int_{\mathbb R} ( \tfrac{m}{m-1}( C - C\tfrac{z^2}{2}  )  )_+^{\frac{1}{m-1}} C^{\frac 1 2} dz \\
#     &= C^{\frac 1 {m-1} + \frac 1 2} \int_{\mathbb R} B_C(z,1) dz.
# \end{aligned}
# $$ 
# The mass of $B_C(\cdot,1)$ can be easily computed

struct Barenblatt{T}
    m::T
    N::Int64
    C::T
    α::T
    β::T
    k::T
end
function Barenblatt(; m, dimension, mass)
    N = dimension
    β = 1 / (2 + N * (m - 1))
    α = β * N
    k = β * (m - 1) / (2 * m)
    basicBarenblatt = Barenblatt{Float64}(m, N, 1.0, α, β, k)
    surfaceArea = 2 * π^(N / 2) / gamma(N / 2)
    integrand(r) = N * surfaceArea * evaluate_radial(basicBarenblatt, r, 1) * r^(N - 1)
    D = quadgk(integrand, 0, Inf, rtol=1e-10)[1]
    C = (mass / D)^(1 / (1 / (m - 1) + 1 / 2))
    return Barenblatt(m, N, C, α, β, k)
end

function evaluate_radial(B::Barenblatt, r::Number, t)
    return t^(-B.α) * max(B.C - B.k * r^2 / t^(2 * B.β), 0.0)^(1 / (B.m - 1))
end

function evaluate(B::Barenblatt, x::Vector, t)
    return evaluate_radial(B, norm(x), t)
end

