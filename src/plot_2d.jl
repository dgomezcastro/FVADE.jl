using Plots, LaTeXStrings

struct UniformMeshADEPlotData2d <: UniformMeshADEPlotData
    i1s
    i2s
    xs
    ys
end

function initialize_plotmesh_2d!(mesh::UniformMeshADE)
    d = dimension(mesh)
    if d ≠ 2
        @error "This function is for d=2"
    end
    imins = [minimum(i[k] for i in mesh.Ih) for k = 1:d]
    imaxs = [maximum(i[k] for i in mesh.Ih) for k = 1:d]

    i1s = imins[1]:imaxs[1]
    i2s = imins[2]:imaxs[2]
    xs = mesh.h[1] .* i1s
    ys = mesh.h[2] .* i2s
    mesh.plotting_object = UniformMeshADEPlotData2d(i1s, i2s, xs, ys)
end

function vector_to_matrix(ρ, mesh::UniformMeshADE)
    d = dimension(mesh)
    if d ≠ 2
        @error "This function is for d=2"
    end

    if isnothing(mesh.plotting_object)
        initialize_plotmesh_2d!(mesh)
    end

    i1s = mesh.plotting_object.i1s
    i2s = mesh.plotting_object.i2s
    ρ_matrix = Matrix(undef, length(i1s), length(i2s))
    ρ_matrix .= Float64(NaN)
    for (p1, i1) in enumerate(i1s)
        for (p2, i2) in enumerate(i2s)
            q = findfirst(==([i1, i2]), mesh.Ih)
            if !isnothing(q)
                ρ_matrix[p1, p2] = ρ[q]
            end
        end
    end
    return ρ_matrix
end

function plot_2d(ρ, mesh; color=nothing, alpha=1.0)
    d = dimension(mesh)
    if d ≠ 2
        @error "This function is for d=2"
    end

    ρ_matrix = vector_to_matrix(ρ, mesh)

    if isnothing(color)
        p = surface(
            mesh.plotting_object.xs,
            mesh.plotting_object.ys,
            ρ_matrix',
            alpha=alpha
        )
    else
        p = surface(
            mesh.plotting_object.xs,
            mesh.plotting_object.ys,
            ρ_matrix',
            c=cgrad([color, color]),
            colorbar=false,
            alpha=alpha
        )
    end

    xlabel!(p, L"x")
    ylabel!(p, L"y")

    return p
end