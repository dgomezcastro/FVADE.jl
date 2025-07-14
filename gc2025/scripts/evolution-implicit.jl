using ProgressMeter
import FVADE
using Plots, LaTeXStrings
using LinearAlgebra

# simulation_name = "driftdiffusion-saturation-square"
# include("cases/$simulation_name.jl")

# ENV["JULIA_DEBUG"] = all
function generate_figure_evolution()

    println("Meshing")
    mesh = FVADE.UniformMeshADE(
        is_in_Omega=is_in_Omega,
        h=h,
        mesh_limits=limits
    )

    println("size Ih = ", length(mesh.Ih))


    ρ = FVADE.initialize_ρ(ρ0, mesh)
    N = ceil(Int64, T / τ) + 1

    println("Solving")
    p = Progress(N)
    anim = @animate for n in 1:N
        FVADE.plot_2d(ρ, mesh)
        if !isnothing(zlimit)
            zlims!(0, zlimit)
        end
        title!("t=$(round((n-1)*τ,digits=3))")

        ρ = FVADE.iterate(
            ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20
        )

        next!(p)
    end
    mp4(anim, "figures/$simulation_name.mp4")

    theplot = FVADE.plot_2d(ρ, mesh)
    title!(plottitle * latexstring(", t=$T"))
    if !isnothing(zlimit)
        zlims!(0, zlimit)
    end
    savefig("figures/$simulation_name.pdf")
end