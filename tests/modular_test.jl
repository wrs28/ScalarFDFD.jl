@startup :orange
push!(LOAD_PATH,"/Users/wrs/wrs julia/ScalarFDFD/src")
using Revise
##

# using SimulationDefinition
# using DifferentialOperators
# using HelmholtzEigen
using ConstructionTools
# using ConstructionTools
# include("../src/helmholtz_overloading.jl")
# include("../src/simulation_plot.jl")
# include("../src/defaults.jl")
##



# L1 = laplacian(11, .1, RobinBoundaryCondition([:d,:n]))
using Arpack

N = [201, 201]
dr = 1/N[1]
dϕ = 2π/N[2]
bc1 = RobinBoundaryCondition(α=[1,1],β=[0,0],p=[0,0],q=[0,0],g=[0,0])
bc2 = PeriodicBoundaryCondition(N,2)
∇², S = laplacian(N, [dr,dϕ], (bc1,bc2); coordinate_system=:polar)
e,v = eigs(∇²; sigma=-10^2)
w = reshape(v[:,5],N[1],N[2])
heatmap(abs.(w)')

# lattice = BravaisLattice(a=1,b=1,β=π/3)
##


##
using Arpack
using Interpolations

e, v = eigs(L2; sigma=3)
w = reshape(v[:,2],N[1],N[2])

using Plots
r = LinRange(0,dx*N[1],N[1])
ϕ = LinRange(0,2π,N[2])

x = r.*cos.(ϕ')
y = r.*sin.(ϕ')

plot(x, y, abs2.(w)', seriestype=:path)

##
@code_wartype PeriodicBoundaryCondition([11,12],1)

##
using Interpolations
W = hcat(w[:,end],w,w[:,1])
r = dr*(1/2 .+ (0:N[1]-1))
ϕ = dϕ*(1/2 .+ (-1:N[2]))
itp = extrapolate(scale(interpolate(spdiagm(0=>1 ./sqrt.(r))*W*conj(W[1,1])/abs(W[1,1]),BSpline(Cubic(Line(OnGrid())))),r,ϕ),NaN*complex(1,1))
x = LinRange(-1,1,maximum(N))
y = LinRange(-1,1,maximum(N))
R = hypot.(x,y')
Φ = mod.(atan.(y',x),2π)
heatmap(x,y,real.(itp.(R,Φ)'))

##
r = dr*(1/2 .+ (0:N[1]-1))
a = w*conj(w[1,1])/abs(w[1,1])
plot(r,real(a[:,1]), xscale=:log10, yscale=:log10)


hcat([0.],)
