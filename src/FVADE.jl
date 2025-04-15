module FVADE

struct ADEProblem
    U::Function
    V::Function
    K::Function
end

include("mesh.jl")

include("equations.jl")

end # module FVADE
