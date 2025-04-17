"""
Mesh for ADE
- h > 0 so that x_i = i * h
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

abstract type MeshADEPlotData end

mutable struct MeshADE
    const h::Float64
    const Ih::Vector{Vector{Int64}}
    const neighbours_plus::Matrix{Union{Int64,Nothing}}
    const neighbours_minus::Matrix{Union{Int64,Nothing}}
    VV::Union{Vector,Nothing}
    KK::Union{Matrix,Nothing}
    plotting_object::Union{MeshADEPlotData,Nothing}
end
function MeshADE(; problem::ADEProblem, is_in_Omega, h, mesh_limits)
    Ih = generate_Ih(is_in_Omega, h, mesh_limits)
    d = length(mesh_limits)
    neighbours_plus = find_neighbours_plus(Ih, d)
    neighbours_minus = find_neighbours_minus(Ih, d)
    VV, KK = initialize_VV_WW(h, Ih, problem.V, problem.K)
    return MeshADE(h, Ih, neighbours_plus, neighbours_minus, VV, KK, nothing)
end
function dimension(mesh::MeshADE)
    length(mesh.Ih[1])
end
function x(i::Vector{Int64}, h::Number)::Vector{Float64}
    return i * h
end

function generate_Ih(is_in_Omega::Function, h::Real, mesh_limits::Vector)::Vector{Vector{Int64}}
    d = length(mesh_limits)
    cube_zip = floor(Int64, mesh_limits[1][1] / h):ceil(Int64, mesh_limits[1][2] / h)
    for k in 2:d
        cube_zip = Iterators.product(cube_zip, floor(Int64, mesh_limits[k][1] / h):ceil(Int64, mesh_limits[k][2] / h))
    end
    Ih = Vector{Vector{Int64}}()
    for i in cube_zip
        if typeof(i) <: Number
            i_vector = [i]
        else
            i_vector = collect(i)
        end
        x_i = x(i_vector, h)
        if is_in_Omega(x_i)
            append!(Ih, [i_vector])
        end
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

function initialize_VV_WW(h, Ih, V, K)
    Nh = length(Ih)
    d = length(Ih[1])
    if isnothing(V)
        VV = nothing
    else
        VV = [V(x(i, h)) for i in Ih]
    end
    if isnothing(K)
        KK = nothing
    else
        KK = Matrix{Float64}(undef, Nh, Nh)
        for (p, i) in enumerate(Ih)
            for (q, j) in enumerate(Ih)
                KK[p, q] = h^d * K(x(i, h), x(j, h))
            end
        end
    end
    return VV, KK
end