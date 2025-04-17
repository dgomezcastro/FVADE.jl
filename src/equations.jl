using NonlinearSolve

function ξ(ρ::Vector{T}, problem::ADEProblem, mesh::MeshADE) where {T<:Number}
    ξ = zeros(T, size(ρ))
    if !(isnothing(problem.Uprime))
        ξ += problem.Uprime.(ρ)
    end
    if !(isnothing(problem.V))
        ξ += mesh.VV
    end
    if !(isnothing(problem.K))
        ξ += mesh.KK * ρ
    end
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
function F(ρ::Vector{T}, v, problem::ADEProblem, mesh::MeshADE) where {T<:Number}
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    F = Matrix{T}(undef, Nh, d)
    mobupwind(a, b) =
        max(
            problem.mobup(a) * problem.mobdown(b),
            0.0
        )
    for p in eachindex(mesh.Ih)
        for k = 1:d
            q = mesh.neighbours_plus[p, k] #i_q = i_p + e_k
            if isnothing(q)
                F[p, k] = 0.0
            elseif v[p, k] ≥ 0
                # Inside implicit solver this values ρ < 0 is possible
                F[p, k] = mobupwind(ρ[p], ρ[q]) * v[p, k]
            elseif v[p, k] < 0
                F[p, k] = mobupwind(ρ[q], ρ[p]) * v[p, k]
            end
        end
    end
    return F
end

"""
L(ρ)[p] = τ * ∑_k (F_{i_p+1/2e_k} - F_{i_p-1/2*e_k})
"""
function L(F::Matrix{T}, mesh, τ) where {T<:Number}
    dρ = zeros(T, length(mesh.Ih))
    for p in eachindex(mesh.Ih)
        for k = 1:dimension(mesh)
            dρ[p] += F[p, k]
            q_minus = mesh.neighbours_minus[p, k] #i_{q_-} = i_p - e_k
            if !(isnothing(q_minus))
                dρ[p] -= F[q_minus, k]
            end
        end
    end
    return (τ / mesh.h) * dρ
end

function implicit_problem(ρ, ρ_prev, F, problem, mesh, τ)
    ξ_ρ = ξ(ρ, problem, mesh)
    v_ρ = v(ξ_ρ, mesh)
    F_ρ = F(ρ, v_ρ, problem, mesh)
    L_ρ = L(F_ρ, mesh, τ)
    return ρ + L_ρ - ρ_prev
end

function iterate(ρ_prev, problem::ADEProblem, mesh::MeshADE, τ::Number; abs_tol=1e-3, max_iters=100)
    G(ρ) = implicit_problem(ρ, ρ_prev, F, problem, mesh, τ)
    ρ_next = Newton(G, ρ_prev; abs_tol=abs_tol, max_iters=max_iters)
    # Newton solve may break positivity and mass conservation
    ρ_next = max.(ρ_next, 0.0)
    if sum(ρ_prev) > 0
        ρ_next = ρ_next / sum(ρ_next) * sum(ρ_prev)
    end
    return ρ_next

    # f(ρ_next, _) = implicit_problem(ρ_next, ρ, F, problem, mesh, τ)
    # prob = NonlinearProblem(f, ρ, [])
    # sol = solve(prob, NewtonRaphson(), abstol=1e-3, maxiters=round(Int64, 1e2))
    # return sol.u
end

function iterate_explicit(ρ_prev, problem::ADEProblem, mesh::MeshADE, τ::Number; abs_tol=1e-3, max_iters=100)
    ξ_ρ = ξ(ρ_prev, problem, mesh)
    v_ρ = v(ξ_ρ, mesh)
    F_ρ = F(ρ_prev, v_ρ, problem, mesh)
    L_ρ = L(F_ρ, mesh, τ)
    ρ_next = ρ_prev - L_ρ
    ρ_next = max.(ρ_next, 0.0)
    if sum(ρ_prev) > 0
        ρ_next = ρ_next / sum(ρ_next) * sum(ρ_prev)
    end
    return ρ_next
end