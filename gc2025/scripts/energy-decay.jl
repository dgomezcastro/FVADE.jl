import FVADE
using ProgressMeter
using Plots, LaTeXStrings
using LinearAlgebra
include("paper-figures.jl")

# experiment = "aggregation-d2"
# include("cases/$experiment.jl")

function generate_figure_energy_decay()
    println("Meshing")
    mesh = FVADE.UniformMeshADE(
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )
    println("size Ih = ", length(mesh.Ih))

    # ENV["JULIA_DEBUG"] = all

    ρ = FVADE.initialize_ρ(ρ0, mesh)

    println("Solving")
    N = ceil(Int64, T / τ)
    p = Progress(N)

    free_energies = Vector{Float64}(undef, N)
    for n in 1:N
        free_energies[n] = FVADE.free_energy(ρ, problem, mesh)

        ρ = FVADE.iterate(
            ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20, solver=paper_solver
        )

        next!(p)
    end

    plot(τ * [0:(N-1);], free_energies,
        label="",
        xlabel=L"t",
        ylabel=L"\mathcal{F}[\rho_t]",
        linewidth=2,
        left_margin=5mm) # so that the label is not cut
    savefig("figures/energydecay-$experiment.pdf")
    save_caption("figures/energydecay-$experiment.tex", plottitle * latexstring(", T=$T"))

    if length(mesh.h) == 1
        theplot = FVADE.plot_1d(ρ, mesh)
        ylabel!(L"\rho_t")
        xlabel!(L"x")
        savefig("figures/energydecay-$experiment-steadystate.pdf")
    elseif length(mesh.h) == 2
        theplot = FVADE.plot_2d(ρ, mesh)
        zlabel!(L"\rho_t")
        label_horizontal_axes!(theplot, mesh.plotting_object.xs, mesh.plotting_object.ys)
        plot!(right_margin=6mm, # so that the colorbar labels are not cut
            bottom_margin=8mm) # and the labels of the horizontal axes
        savefig("figures/energydecay-$experiment-steadystate.pdf")
    end
    save_caption("figures/energydecay-$experiment-steadystate.tex", plottitle * latexstring(", t=$T"))
end
