using Plots, LaTeXStrings

function plot_1d(ρ, mesh; color=nothing, alpha=1.0, label="")
    d = dimension(mesh)
    if d ≠ 1
        @error "This function is for d=1"
    end

    x = [FVADE.x(i, mesh.h)[1] for i in mesh.Ih]

    if isnothing(color)
        p = plot(
            x, ρ,
            label=label,
            alpha=alpha
        )
    else
        p = plot(
            x, ρ,
            label=label,
            color=color,
            alpha=alpha
        )
    end

    return p
end

function plot_1d!(ρ, mesh; color=nothing, alpha=1.0, label="")
    d = dimension(mesh)
    if d ≠ 1
        @error "This function is for d=1"
    end

    x = [FVADE.x(i, mesh.h)[1] for i in mesh.Ih]

    if isnothing(color)
        p = plot!(
            x, ρ,
            label=label,
            alpha=alpha
        )
    else
        p = plot!(
            x, ρ,
            label=label,
            color=color,
            alpha=alpha
        )
    end

    return p
end
