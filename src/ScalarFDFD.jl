module ScalarFDFD

using Arpack,
Distributed,
Formatting,
Interpolations,
IterTools,
LinearAlgebra,
LaTeXStrings,
NLsolve,
Optim,
ProgressMeter,
Random,
RecipesBase,
Serialization,
SparseArrays,
SpecialFunctions,
Statistics

# import Optim: optimize, BFGS
# import UnicodePlots: lineplot, lineplot!, BlockCanvas

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

# Eigenvalue Solvers, all checked in 1d (except β)
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
export band_structure

# Waveguide Dispersions
export build_dispersions!,
remove_dispersions!,
waveguide_dispersion

# Plotting and Animation # 1d checked
export wave

# Parallel
export eig_kp,
band_structure_p,
build_dispersions_p!,
waveguide_dispersion_p,
smatrix_p

# Front-end Simulation Construction
export Regions,
build_domain,
add_regions,
add_region,
add_circle_region,
add_circle_domain,
add_ellipse,
add_square,
add_rectangle_region,
add_rectangle_domain,
add_parallelogram

export build_pc_domain,
defect_domain,
add_defect,
line_defect_domain,
add_circle_to_pc

export add_planar_waveguide,
add_planar_waveguides,
add_horizontal_planar_waveguides,
add_vertical_planar_waveguides, #not checked, don't care for now
add_halfspace_waveguide,
add_halfspace_waveguides,
add_pc_waveguide #not checked, don't care for now


# Simulation Definitions
include("defaults.jl")
include("definitions.jl")
include("expanded definitions.jl")
include("iss.jl")
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
include("dispersions.jl")
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
include("construction_framework.jl")
include("shapes.jl")
include("shape_constructors.jl")
include("pc.jl")
include("waveguides.jl")

end
