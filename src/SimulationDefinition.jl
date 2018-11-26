module SimulationDefinition

using Interpolations

export Bravais,
    Domain,
    System,
    Discretization,
    Boundary,
    Channels,
    Scattering,
    TwoLevelSystem,
    Simulation

    include("defaults.jl")
    include("simulation structures.jl")
    include("simulation construction.jl")
    include("simulation overloading.jl")
    include("bravais coordinates.jl")

    include("iss.jl")

end # module
