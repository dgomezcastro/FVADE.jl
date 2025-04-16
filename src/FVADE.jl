module FVADE

using ForwardDiff

struct ADEProblem
    U::Function
    Uprime::Function
    V::Function
    K::Union{Function,Nothing}
end
function ADEProblem(U, V, K; Uprime=s -> ForwardDiff.derivative(U, s))
    return ADEProblem(U, Uprime, V, K)
end

include("mesh.jl")

include("solver.jl")
include("equations.jl")

include("plot_2d.jl")

end # module FVADE
