using Plots, Plots.Measures, LaTeXStrings

# The figures are included in the paper with width 0.45\textwidth, i.e., scaled down to about a
# quarter of their size, so the labels need large fonts to be readable
default(linewidth=2, guidefontsize=20, tickfontsize=20, legendfontsize=18, colorbar_tickfontsize=20)

# The numbers in the axes are written in LaTeX, like the labels. GR ignores this for the colorbar, so we
# use the font of LaTeX for the rest of the text, and a smaller size, since GR draws LaTeX smaller
default(formatter=:latex, fontfamily="Computer Modern", colorbar_tickfontsize=16)

mkpath("figures")

# Solver for FVADE.iterate. The figures in the paper use :NewtonRaphson, :SimpleDFSane is much faster
paper_solver = :NewtonRaphson

# Writes the parameters in `plottitle` as a LaTeX snippet, to be \input in a (sub)caption
function save_caption(filename, plottitle)
    write(filename, replace(plottitle, "\$\n\$" => ", ", "\$\$" => "", "τ" => "\\tau") * "\n")
end

# Limits of a logarithmic y axis that leave room for a legend above the curves ys
function legend_room(ys...)
    lo, hi = minimum(minimum.(ys)), maximum(maximum.(ys))
    return (lo / 1.2, hi * (hi / lo)^0.4)
end

# GR draws the tick labels of 3D plots very close to the axes, so we pad those of the vertical axis with spaces
padded_ticks(ticks) = (ticks, [latexstring("$(t)\\;\\;") for t in ticks])

# At most three ticks on each horizontal axis, so that the labels of the 3D plots do not overlap
function sparse_ticks(xs)
    ticks = (2*ceil(Int, minimum(xs) / 2)):2:(2*floor(Int, maximum(xs) / 2))
    return length(ticks) > 3 ? ticks[1:2:end] : ticks
end

# GR does not allow to move the tick labels of the horizontal axes of 3D plots away from the axes, so
# we write them, and the names of the axes, as annotations at a distance proportional to the box
function label_horizontal_axes!(p, xs, ys)
    xticks, yticks = sparse_ticks(xs), sparse_ticks(ys)
    dx, dy = 0.14 * (maximum(xs) - minimum(xs)), 0.14 * (maximum(ys) - minimum(ys))
    plot!(p, xlims=extrema(xs), ylims=extrema(ys), xlabel="", ylabel="",
        xticks=(xticks, fill("", length(xticks))), yticks=(yticks, fill("", length(yticks))))
    number(t) = text(latexstring(string(t)), Plots.default(:tickfontsize), :top)
    for t in xticks
        annotate!(p, t, minimum(ys) - dy, 0.0, number(t))
    end
    for t in yticks
        annotate!(p, maximum(xs) + dx, t, 0.0, number(t))
    end
    annotate!(p, 0.0, minimum(ys) - 3dy, 0.0, text(L"x", Plots.default(:guidefontsize), :top))
    annotate!(p, maximum(xs) + 3dx, 0.0, 0.0, text(L"y", Plots.default(:guidefontsize), :top))
    return p
end
