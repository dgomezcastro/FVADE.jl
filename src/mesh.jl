using LinearAlgebra

abstract type UniformMeshADEPlotData end
abstract type ProblemCoefficients end

"""
Mesh for ADE
    - h vector so that x_i = i * h[i]
    - Ih = [i ∈ Z^d such that x_i ∈ Ω]
    - neighbours_plus[p,k] = q 
    - if q∈N then x_{i_q} = x_{i_p + e_k} ∈ Ω
    - q=nothing then x_{i_p + e_k} ∉ Ω
    - neighbours_minus[p,k] = q 
    - if q∈N then x_{i_q} = x_{i_p - e_k} ∈ Ω
    - q=nothing then x_{i_p - e_k} ∉ Ω
    
    - VV[p] = V(x[p])
    - KK[p,q] = K(x[p],x[q])
"""
mutable struct UniformMeshADE
    const h::Vector{Float64}
    const Ih::Vector{Vector{Int64}}
    Ihplus::Union{Vector{Vector{Int64}},Nothing} # Additional space to add more points to the mesh
    const neighbours_plus::Matrix{Union{Int64,Nothing}}
    const neighbours_minus::Matrix{Union{Int64,Nothing}}
    coefficients::Union{ProblemCoefficients,Nothing}
    plotting_object::Union{UniformMeshADEPlotData,Nothing}
end
function UniformMeshADE(; is_in_Omega, h::Union{Real,Vector}, mesh_limits)
    d = length(mesh_limits)
    if typeof(h) <: Number
        h = h * ones(d)
    end
    @debug typeof(h)

    Ih = generate_Ih(is_in_Omega, h, mesh_limits)
    d = length(mesh_limits)
    neighbours_plus = find_neighbours_plus(Ih, d)
    neighbours_minus = find_neighbours_minus(Ih, d)
    return UniformMeshADE(h, Ih, nothing, neighbours_plus, neighbours_minus, nothing, nothing)
end

function dimension(mesh::UniformMeshADE)
    length(mesh.h)
end
function cubevolume(h::Vector)
    return prod(h for h in h)
end

function x(i::Vector{Int64}, h::Vector)::Vector{Float64}
    return i .* h
end

function x(i::Vector{Int64}, h::Number)::Vector{Float64}
    return i * h
end

function generate_Ih(is_in_Omega::Function, h::Vector, mesh_limits::Vector)::Vector{Vector{Int64}}
    d = length(mesh_limits)
    if d == 1
        range = floor(Int64, mesh_limits[1][1] / h[1]):ceil(Int64, mesh_limits[1][2] / h[1])
        Ih = [[i] for i in range if is_in_Omega(x([i], h))]
    else
        cube_zip = floor(Int64, mesh_limits[1][1] / h[1]):ceil(Int64, mesh_limits[1][2] / h[1])
        for k in 2:d
            cube_zip = Iterators.product(cube_zip, floor(Int64, mesh_limits[k][1] / h[k]):ceil(Int64, mesh_limits[k][2] / h[k]))
        end
        Ih = [collect(i) for i in cube_zip if is_in_Omega(x(collect(i), h))]
    end
    return Ih
end

function find_neighbours_plus(Ih, d)::Matrix{Union{Int64,Nothing}}
    Nh = length(Ih)
    neighbours = Matrix{Union{Int64,Nothing}}(undef, Nh, d)
    for (p, i) in enumerate(Ih)
        for k = 1:d
            ek = zeros(Int64, d)
            ek[k] = 1
            neighbours[p, k] = findfirst(==(i + ek), Ih)
        end
    end
    return neighbours
end


function find_neighbours_minus(Ih, d)::Matrix{Union{Int64,Nothing}}
    Nh = length(Ih)
    neighbours = Matrix{Union{Int64,Nothing}}(undef, Nh, d)
    for (p, i) in enumerate(Ih)
        for k = 1:d
            ek = zeros(Int64, d)
            ek[k] = 1
            neighbours[p, k] = findfirst(==(i - ek), Ih)
        end
    end
    return neighbours
end

function initialize_ρ(ρ0::Function, mesh)::Vector
    return [Float64(ρ0(x(i, mesh.h))) for i in mesh.Ih]
end