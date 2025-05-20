using LinearAlgebra

struct MeshCoefficientsADE <: ProblemCoefficients
    VV::Union{Vector,Nothing}
    KK::Union{Symmetric{Float64,Matrix{Float64}},Nothing}
end

function initialize_VV_WW(h, Ih, V, K)::MeshCoefficientsADE
    if isnothing(V)
        VV = nothing
    else
        VV = [V(x(i, h)) for i in Ih]
    end
    if isnothing(K)
        KK = nothing
    else
        KK = Symmetric([cubevolume(h) * K(x(i, h), x(j, h)) for i in Ih, j in Ih])
    end
    return MeshCoefficientsADE(VV, KK)
end

function initialize!(mesh::UniformMeshADE, problem::ADEProblem)
    mesh.coefficients = initialize_VV_WW(mesh.h, mesh.Ih, problem.V, problem.K)
    return nothing
end

function ξ(ρ_next::Vector{T}, ρ_prev, problem::ADEProblem, mesh::UniformMeshADE) where {T<:Number}
    ξ = zeros(T, size(ρ_next))
    if !(isnothing(problem.Uprime))
        ξ += problem.Uprime.(max.(ρ_next, 0.0))
    end
    if !(isnothing(problem.V))
        ξ += mesh.coefficients.VV
    end
    if !(isnothing(problem.K))
        ξ += mesh.coefficients.KK * (ρ_next + ρ_prev) / 2
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

"""
Iterates the ADE problem with no-flux conditions using implicit Euler
"""
function iterate(ρ_prev, problem::ADEProblem, mesh::UniformMeshADE, τ::Number; abs_tol=1e-3, max_iters=100)
    if isnothing(mesh.coefficients)
        initialize!(mesh, problem)
    end

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

"""
Iterates the ADE problem with no-flux conditions using explicit Euler
"""
function iterate_explicit(ρ_prev, problem::ADEProblem, mesh::UniformMeshADE, τ::Number)
    if isnothing(mesh.coefficients)
        initialize!(mesh, problem)
    end

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
        free_energy += vol * dot(mesh.coefficients.VV, ρ)
    end
    if !isnothing(problem.K)
        free_energy += 0.5 * vol * dot(ρ, (mesh.coefficients.KK) * ρ)
    end
    return free_energy
end