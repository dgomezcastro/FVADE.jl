include("../scripts/energy-decay.jl")

experiment = "aggregation-d1"
include("../scripts/cases/$experiment.jl")
generate_figure_energy_decay()

experiment = "aggregation-d2"
include("../scripts/cases/$experiment.jl")
generate_figure_energy_decay()
