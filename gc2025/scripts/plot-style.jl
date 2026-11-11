using Plots, Plots.Measures

# The figures are included in the paper with width 0.45\textwidth (≈155pt), so the canvas is
# small and the fonts large: once scaled down, titles and ticks are ≈7pt and axis labels ≈8.5pt.
default(
    size=(380, 300),
    titlefontsize=13,
    guidefontsize=16,
    tickfontsize=13,
    legendfontsize=12,
    linewidth=2,
    markersize=4,
    margin=3mm,
    top_margin=6mm, # the titles have three lines
)
