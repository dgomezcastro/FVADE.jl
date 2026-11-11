using ProgressMeter
import FVADE
using Plots, LaTeXStrings
using LinearAlgebra
include("paper-figures.jl")

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

    # Axes of the 3D plots, both in the video and in the figure
    function style_axes!(theplot)
        if !isnothing(zlimit)
            plot!(theplot, zlims=(0, zlimit), zticks=padded_ticks(range(0, zlimit, length=3)))
        end
        label_horizontal_axes!(theplot, mesh.plotting_object.xs, mesh.plotting_object.ys)
        plot!(theplot, right_margin=6mm, # so that the colorbar labels are not cut
            bottom_margin=8mm) # and the labels of the horizontal axes
    end

    println("Solving")
    p = Progress(N)
    anim = @animate for n in 1:N
        style_axes!(FVADE.plot_2d(ρ, mesh))
        plot!(title="t=$(round((n-1)*τ,digits=3))", titlefontsize=20)

        ρ = FVADE.iterate(
            ρ, problem, mesh, τ; abs_tol=1e-8, max_iters=20, solver=paper_solver
        )

        next!(p)
    end
    # mp4(anim, "figures/$simulation_name.mp4")

    style_axes!(FVADE.plot_2d(ρ, mesh))
    savefig("figures/$simulation_name.pdf")
    save_caption("figures/$simulation_name.tex", plottitle * latexstring(", t=$T"))
end