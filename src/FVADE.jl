module FVADE

using NonlinearSolve

struct ADEProblem
    U::Function
    V::Function
    K::Function
end

include("mesh.jl")

function ξ(ρ, problem::ADEProblem, mesh::MeshADE)
    cube_volume = (mesh.h)^(dimension(mesh))
    ξ = [
        problem.U(ρ[p]) +
        mesh.VV[p] +
        cube_volume * sum(mesh.KK[p, q] * ρ[q] for q in eachindex(mesh.Ih))
        for p in eachindex(mesh.Ih)
    ]
    return ξ
end

"""
v[p,k] = v_{i_p + 1/2 * e_k}    
"""
function v(ξ, mesh::MeshADE)
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    v = Matrix{Float64}(undef, Nh, d)
    for p in eachindex(mesh.Ih)
        for k = 1:d
            q = mesh.neighbours_plus[p, k] #i_q = i_p + e_k
            if isnothing(q)
                v[p, k] = 0.0
            else
                v[p, k] = -(ξ[q] - ξ[p]) / (mesh.h)
            end
        end
    end
    return v
end

"""
F[p,k] = F_{i_p + 1/2 * e_k}    
"""
function F(ρ, v, mesh::MeshADE)
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    F = Matrix{Float64}(undef, Nh, d)
    for p in eachindex(mesh.Ih)
        for k = 1:d
            q = mesh.neighbours_plus[p, k] #i_q = i_p + e_k
            if isnothing(q)
                F[p, k] = 0.0
            elseif v[p, k] ≥ 0
                F[p, k] = ρ[p] * v[p, k]
            else
                F[p, k] = ρ[q] * v[p, k]
            end
        end
    end
    return v
end

"""
H(ρ)[p] = ρ_{i_p} + τ * ∑_k (F_{i_p+1/2e_k} - F_{i_p-1/2*e_k})
"""
function H(ρ, F, mesh, τ)
    dρ = zeros(length(mesh.Ih))
    for p in eachindex(mesh.Ih)
        for k = 1:dimension(mesh)
            q_plus = mesh.neighbours_plus[p, k] #i_q = i_p + e_k
            q_minus = mesh.neighbours_minus[p, k] #i_q = i_p - e_k
            if !(isnothing(q_plus))
                dρ[p] += F[p, k]
            end
            if !(isnothing(q_minus))
                dρ[p] -= F[q_minus, k]
            end
        end
    end
    return ρ + (τ / mesh.h) * dρ
end

function iterate(ρ, problem, mesh::MeshADE, τ)
    ξ_at_ρ = ξ(ρ, problem, mesh)
    v_at_ρ = v(ξ_at_ρ, mesh)
    F_at_ρ = F(ρ, v_at_ρ, mesh)

    f(ρ, _) = H(ρ, F_at_ρ, mesh, τ)
    prob = NonlinearProblem(f, ρ, [])
    sol = solve(prob, NewtonRaphson())
    return sol.u
end



end # module FVADE
