module FVADE

using ForwardDiff

struct ADEProblem
    U::Union{Function,Nothing}
    Uprime::Union{Function,Nothing}
    V::Union{Function,Nothing}
    K::Union{Function,Nothing}
    mobup::Function
    mobdown::Function
end
function ADEProblem(;
    U::Union{Function,Nothing}=nothing,
    V::Union{Function,Nothing}=nothing,
    K::Union{Function,Nothing}=nothing,
    mobup::Function=s -> s,
    mobdown::Function=s -> 1,
    Uprime=isnothing(U) ? nothing : (s -> ForwardDiff.derivative(U, s)))
    return ADEProblem(U, Uprime, V, K, mobup, mobdown)
end

include("mesh.jl")

include("solver.jl")
include("equations.jl")

include("Barenblatt.jl")

include("plot_1d.jl")
include("plot_2d.jl")

end # module FVADE
