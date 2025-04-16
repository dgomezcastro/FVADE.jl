using Plots

function plot_2d(ρ, mesh)
    d = dimension(mesh)
    if d ≠ 2
        @error "This function is for d=2"
    end

    imins = [minimum(i[k] for i in mesh.Ih) for k = 1:d]
    imaxs = [maximum(i[k] for i in mesh.Ih) for k = 1:d]

    i1s = imins[1]:imaxs[1]
    xs = mesh.h * i1s
    i2s = imins[2]:imaxs[2]
    ys = mesh.h * i2s

    ρ_matrix = zeros(length(i1s), length(i2s))
    for (p1, i1) in enumerate(i1s)
        for (p2, i2) in enumerate(i2s)
            q = findfirst(==([i1, i2]), mesh.Ih)
            if !isnothing(q)
                ρ_matrix[p1, p2] = ρ[q]
            end
        end
    end

    surface(xs, ys, ρ_matrix)
end