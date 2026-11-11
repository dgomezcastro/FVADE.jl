# Plots figure 4 from the errors saved in figures/*.jld2 by figure4-convergence-halfh.jl
using Plots, LaTeXStrings, Plots.Measures
using JLD2
include("scripts/paper-figures.jl")

function plot_figure4(name)
    data = load("figures/$name.jld2")
    hs, L1errors = data["hs"], data["L1errors"]

    plot(hs, L1errors,
        scale=:log10,
        marker=:circle,
        xlabel=L"h",
        ylabel=L"ε_h^{(2)}",
        label="Simulation",
        size=(700, 400),
        bottommargin=5mm,
        leftmargin=5mm,
        legend=:topleft,
        ylims=legend_room(L1errors, hs)
    )
    plot!(hs, hs, label=L"Linear scaling $\varepsilon_h = h$")

    savefig("figures/$name.pdf")
    if haskey(data, "plottitle")
        save_caption("figures/$name.tex", data["plottitle"])
    else
        @warn "figures/$name.jld2 has no plottitle, so figures/$name.tex is not written"
    end
end

for d in (1, 2)
    if isfile("figures/convergence-hhalf-d$d.jld2")
        plot_figure4("convergence-hhalf-d$d")
    elseif d == 1
        @warn "figures/convergence-hhalf-d$d.jld2 not found, run figure4-convergence-halfh.jl"
    end
end
