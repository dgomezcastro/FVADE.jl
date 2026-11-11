# Plots figure 3 from the errors saved in figures/*.jld2 by figure3-convergence-Barenblatt.jl
using Plots, LaTeXStrings, Plots.Measures
using JLD2
include("scripts/paper-figures.jl")

function plot_figure3(name)
    data = load("figures/$name.jld2")
    hs, L1errors = data["hs"], data["L1errors"]

    # Reference line C h, just above the data
    C = 1.5 * maximum(L1errors ./ hs)

    plot(hs, L1errors,
        scale=:log10,
        marker=:circle,
        xlabel=L"h",
        ylabel=L"\varepsilon^{(1)}_h",
        label="Simulation",
        size=(700, 400),
        bottommargin=5mm,
        leftmargin=5mm,
        legend=:topleft,
        ylims=legend_room(L1errors, C * hs)
    )
    plot!(hs, C * hs, label=L"$\varepsilon_h = C h$")

    savefig("figures/$name.pdf")
    if haskey(data, "plottitle")
        save_caption("figures/$name.tex", data["plottitle"])
    else
        @warn "figures/$name.jld2 has no plottitle, so figures/$name.tex is not written"
    end
end

for d in (1, 2)
    if isfile("figures/convergence-Barenblatt-$(d)d.jld2")
        plot_figure3("convergence-Barenblatt-$(d)d")
    else
        @warn "figures/convergence-Barenblatt-$(d)d.jld2 not found, run figure3-convergence-Barenblatt.jl with d = $d"
    end
end
