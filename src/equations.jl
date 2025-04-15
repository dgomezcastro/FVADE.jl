using NonlinearSolve

function ξ(ρ::Vector{T}, problem::ADEProblem, mesh::MeshADE) where {T<:Number}
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
function v(ξ::Vector{T}, mesh::MeshADE) where {T<:Number}
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    v = Matrix{T}(undef, Nh, d)
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
function F(ρ::Vector{T}, v, mesh::MeshADE) where {T<:Number}
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    F = Matrix{T}(undef, Nh, d)
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
function H(ρ::Vector{T}, F, mesh, τ) where {T<:Number}
    dρ = zeros(T, length(mesh.Ih))
    for p in eachindex(mesh.Ih)
        for k = 1:dimension(mesh)
            q_plus = mesh.neighbours_plus[p, k] #i_{q_+} = i_p + e_k
            q_minus = mesh.neighbours_minus[p, k] #i_{q_-} = i_p - e_k
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

function implicit_problem(ρ_next, ρ, F, problem, mesh, τ)
    ξ_next = ξ(ρ_next, problem, mesh)
    v_next = v(ξ_next, mesh)
    F_next = F(ρ_next, v_next, mesh)
    H_next = H(ρ_next, F_next, mesh, τ)
    return H_next - ρ
end

function iterate(ρ, problem, mesh::MeshADE, τ)
    f(ρ_next, _) = implicit_problem(ρ_next, ρ, F, problem, mesh, τ)
    prob = NonlinearProblem(f, ρ, [])
    sol = solve(prob, NewtonRaphson())
    return sol.u
end