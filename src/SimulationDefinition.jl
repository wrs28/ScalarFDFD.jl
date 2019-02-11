module SimulationDefinition

using Bravais,
CoordinateSystems,
BoundaryConditions,
DielectricFunctions,
Formatting,
Interpolations,
IterTools,
RecipesBase,
Shapes,
Statistics

import BoundaryConditions: reorder, get_dim, get_side, apply_args
import DifferentialOperators: _oc_bls

export Domain,
System,
Discretization,
Boundary,
Channels,
Scattering,
TwoLevelSystem,
Simulation,
which_domain

include("defaults.jl")
include("simulation_types.jl")
include("simulation_construction.jl")

include("simulation_overloading.jl")
include("coordinate_overloading.jl")
include("bravais_overloading.jl")

include("iss.jl")

include("simulation_pretty_printing.jl")
include("simulation_plot_defaults.jl")
include("simulation_plot.jl")

end # module
