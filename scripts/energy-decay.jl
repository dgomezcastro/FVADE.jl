import FVADE
using ProgressMeter
using Plots, LaTeXStrings
using LinearAlgebra

# experiment = "aggregation-d2"
# include("cases/$experiment.jl")

function generate_figure_energy_decay()
    println("Meshing")
    mesh = FVADE.UniformMeshADE(
        problem=problem,
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
            ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
        )

        next!(p)
    end

    plot(τ * [0:(N-1);], free_energies,
        title=plottitle * latexstring("T=$T"),
        label="",
        xlabel=L"t",
        ylabel=L"\mathcal{F}[\rho_t]",
        linewidth=2)
    savefig("figures/energydecay-$experiment.pdf")
end
