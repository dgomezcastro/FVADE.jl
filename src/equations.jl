using LinearAlgebra

function ξ(ρ_next::Vector{T}, ρ_prev, problem::ADEProblem, mesh::UniformMeshADE) where {T<:Number}
    ξ = zeros(T, size(ρ_next))
    if !(isnothing(problem.Uprime))
        ξ += problem.Uprime.(max.(ρ_next, 0.0))
    end
    if !(isnothing(problem.V))
        ξ += mesh.VV
    end
    if !(isnothing(problem.K))
        ξ += mesh.KK * (ρ_next + ρ_prev) / 2
    end
    return ξ
end

"""
v[p,k] = v_{i_p + 1/2 * e_k}    
"""
function v(ξ::Vector{T}, mesh::UniformMeshADE) where {T<:Number}
    Nh = length(mesh.Ih)
    d = dimension(mesh)
    v = Matrix{T}(undef, Nh, d)
    for p in eachindex(mesh.Ih)
        for k = 1:d
            q = mesh.neighbours_plus[p, k] #i_q = i_p + e_k
            if isnothing(q)
                v[p, k] = 0.0
            else
                v[p, k] = -(ξ[q] - ξ[p]) / (mesh.h[k])
            end
        end
    end
    return v
end

"""
F[p,k] = F_{i_p + 1/2 * e_k}    
"""
function F(ρ::Vector{T}, v, problem::ADEProblem, mesh::UniformMeshADE) where {T<:Number}
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
            dρ[p] += F[p, k] / mesh.h[k]
            q_minus = mesh.neighbours_minus[p, k] #i_{q_-} = i_p - e_k
            if !(isnothing(q_minus))
                dρ[p] -= F[q_minus, k] / mesh.h[k]
            end
        end
    end
    return τ * dρ
end

function implicit_problem(ρ_next, ρ_prev, problem, mesh, τ)
    ξ_ρ = ξ(ρ_next, ρ_prev, problem, mesh)
    v_ρ = v(ξ_ρ, mesh)
    F_ρ = F(ρ_next, v_ρ, problem, mesh)
    L_ρ = L(F_ρ, mesh, τ)
    return ρ_next + L_ρ - ρ_prev
end

function iterate(ρ_prev, problem::ADEProblem, mesh::UniformMeshADE, τ::Number; abs_tol=1e-3, max_iters=100)
    G(ρ) = implicit_problem(ρ, ρ_prev, problem, mesh, τ)
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

function iterate_explicit(ρ_prev, problem::ADEProblem, mesh::UniformMeshADE, τ::Number)
    ξ_ρ = ξ(ρ_prev, ρ_prev, problem, mesh)
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

function free_energy(ρ, problem, mesh)
    free_energy = 0.0
    vol = cubevolume(mesh.h)
    if !isnothing(problem.U)
        free_energy += vol * sum(problem.U.(ρ))
    end
    if !isnothing(problem.V)
        free_energy += vol * dot(mesh.VV, ρ)
    end
    if !isnothing(problem.K)
        free_energy += 0.5 * vol * dot(ρ, (mesh.KK) * ρ)
    end
    return free_energy
end