include("scripts/evolution-implicit.jl")

zlimit = nothing
simulation_name = "driftdiffusion-saturation-peanut"
include("scripts/cases/$simulation_name.jl")
generate_figure_evolution()

zlimit = nothing
simulation_name = "driftdiffusion-saturation-square"
include("scripts/cases/$simulation_name.jl")
generate_figure_evolution()