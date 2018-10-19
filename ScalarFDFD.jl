module ScalarFDFD

using Arpack,
Distributed,
Formatting,
Interpolations,
IterTools,
LinearAlgebra,
LaTeXStrings,
NLsolve,
ProgressMeter,
Optim,
Random,
RecipesBase,
Serialization,
SparseArrays,
SpecialFunctions,
Statistics

import UnicodePlots: lineplot, lineplot!, BlockCanvas

# Simulation Definitions
export Bravais,
Domain, # incomplete documentation
System,
Discretization,
Boundary,
Channels,
Scattering,
TwoLevelSystem,
Simulation,
extract_waveguide_simulation

# Eigenvalue Solvers
export eig_k,
eig_cf,
eig_β, # did not check, want to revisit this solution
quadrature

# Scattering Solvers
export scattering,
surface_flux, # fails unitarity
bulk_absorption, # fails unitarity
smatrix #fails unitarity

# Band Structure
export band_structure,
build_dispersion!,
remove_dispersion!,
waveguide_dispersion

# Plotting and Animation
export wave

# Parallel
export eig_kp,
band_structure_p,
build_dispersion_p!,
waveguide_dispersion_p,
smatrix_p

# Front-end Simulation Construction
export Regions,
build_domain,
add_regions,
add_region,
add_circle,
add_ellipse,
add_square,
add_rectangle,
add_parallelogram,
build_pc,
add_defect,
add_line_defect,
add_circle,
planar_waveguide, #not checked, don't care for now
add_planar_waveguide, #not checked, don't care for now
defect_waveguide, #not checked, don't care for now
add_defect_waveguide #not checked, don't care for now


# Simulation Definitions
include("defaults.jl")
include("definitions.jl")
include("expanded definitions.jl")
include("bravais.jl")
include("simulation construction.jl")
include("differential_operators.jl")
include("SALT.jl")

# Eigenvalue Solvers
include("eigensolver wrappers.jl")
include("linear eigensolvers.jl")
include("nonlinear eigensolvers.jl")
include("quadrature.jl")

# Scattering Solvers
include("scattering.jl")
include("sources.jl")
include("analysis.jl")
include("smatrix.jl")

# Band Structure
include("band_structure.jl")

# Plotting and Animation
include("pretty_printing.jl")
include("plot.jl")

# Parallel
include("parallel.jl")

# Front-end Simulation Construction
include("shapes.jl")
include("construction_framework.jl")
include("constructors.jl")
include("pc.jl")
include("waveguides.jl")

end
