module StochasticMortalityModels



include("requirements.jl")

include("hmd/hmd.jl")
include("utilities/utilities.jl")
include("models/models.jl")


export greet
greet() = println("Hello from StochasticMortalityModels V2")

end
