# TODO: Test arbitrary profile waveguide, for now have only tested simple waveguide
module Waveguides

export add_halfspace_waveguide,
add_halfspace_waveguides,
add_planar_waveguide,
add_planar_waveguides,
add_pc_waveguide,
add_pc_waveguides

using ...Bravais,
...BoundaryConditions,
...Shapes,
...SimulationDefinition,
..ConstructionToolsBase,
..PhotonicCrystal

include("waveguides_planar.jl")
include("waveguides_halfspace.jl")
include("waveguides_pc.jl")

include("iss.jl")

end # module
