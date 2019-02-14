module ScalarFDFD

include("Bravais.jl"); using .Bravais
export BravaisLattice

include("CoordinateSystems.jl"); using .CoordinateSystems
export Polar,
Cartesian

include("BoundaryConditions.jl"); using .BoundaryConditions
export PML,
cPML,
noBL,
DirichletBC,
NeumannBC,
RobinBC,
PeriodicBC,
MatchedBC

include("DielectricFunctions.jl"); using .DielectricFunctions

include("Shapes.jl"); using .Shapes
export Circle,
Ellipse,
Square,
Rectangle,
Parallelogram,
Universe

include("DifferentialOperators.jl"); using .DifferentialOperators

include("SimulationDefinition.jl"); using .SimulationDefinition
export Domain,
System,
Discretization,
Boundary,
Channels,
Scattering,
TwoLevelSystem,
Simulation

include("ConstructionTools.jl"); using .ConstructionTools
export Subdomains,
add_subdomain,
build_domain,
build_pc_domain,
line_defect,
site_defect,
add_halfspace_waveguide,
add_halfspace_waveguides,
add_planar_waveguide,
add_planar_waveguides,
add_pc_waveguide,
add_pc_waveguides

include("HelmholtzEigen.jl"); using .HelmholtzEigen
export eig_kl,
eig_cf,
eig_knl

end # module
