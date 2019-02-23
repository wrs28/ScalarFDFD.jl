using Plots
if nprocs()==1
      addprocs(7)
      @everywhere using Pkg
      @everywhere Pkg.activate("/Users/wrs/wrs julia/ScalarFDFD")
      @everywhere push!(LOAD_PATH,"/Users/wrs/wrs julia/ScalarFDFD/src")
end
using Revise
using ScalarFDFD
# using BenchmarkTools


##

pc = build_pc_domain(BravaisLattice(a=1, b=1, β=π/3))
pc = Circle(.3, :center, Dict(:n1=>3), pc; x0=.1, y0=.1)
pc = Circle(.1, :center, Dict(:n1=>2), pc)
sys = System(pc)

sys = add_pc_waveguide(sys; width=1, direction=:vertical, x0=0, y0=0)
sys = add_pc_waveguide(sys; width=1, direction=:horizontal, x0=0, y0=0)

bnd = Boundary(∂Ω=(-3,3,-3,3), bc=PeriodicBC)

dis = Discretization(.025)

sim = Simulation(sys, bnd, dis)

plot(sim)
##

bnd = Boundary(∂Ω = (0,π,0,1.1π),bc = DirichletBC)
sys = System(Universe(1))
dis = Discretization(.05)
sim = Simulation(bnd,sys,dis)
##

bnd = Boundary(bl=PML, depth=.1, r=1)
circ = Circle(.5,0,0,Dict(:n1=>3))
sys = System(circ,Universe(1))
dis = Discretization((.01,1*π/180),Polar)

sim = Simulation(bnd,sys,dis)

plot(sim)
##

bnd = Boundary(bc=MatchedBC, r=1.5, outgoing_qns=-5:-5, incoming_qns=-4:5)
circ1 = Circle(.5,0,0,Dict(:n1=>3))
circ2 = Circle(.2,.1,.1,Dict(:n1=>1))
sys = System(circ2,circ1,Universe(1))
dis = Discretization((.02,2*π/180),Polar)

sim = Simulation(bnd,sys,dis)

plot(sim)

##

k, ψ = eig_knl(sim, 3.15; radius=(.75,.75), nk=10, quad_n=50, disp_opt=true)

plot(sim, ψ)
##
##
bc = (DirichletBC{1,2}(),DirichletBC{1,1}(),PeriodicBC{2,2}(),PeriodicBC{2,1}())
lattice = BravaisLattice(b=2π)
bl = (PML{1,2}(),)
N = (600,600)
Δ = (0, 3, 0, 2π)

k = 4.28
σ = k^2
@time (∇², fs, s), defs = laplacian(N,Polar,Δ,bc,bl,k,0,0,1; depth=.5, lattice=lattice);
A = sum(∇²)
##
m = 4
ε = ones(ComplexF64,N[1],N[2])
# ε = ones(ComplexF64,N[1])
ε[1:N[1]÷2,:] .= 9.0
ε = spdiagm(0=>ε[:])
# ε = kron(sparse(I,N[2],N[2]),ε)
r = LinRange(defs[1].Δ[1][1]+defs[1].dx[1]/2,defs[1].Δ[1][2]-defs[1].dx[1]/2,N[1])
# decomp, history = partialschur(shift_and_invert(A-s*m^2*spdiagm(0=>1 ./r.^2),-s*ε,σ),nev=2);
decomp, history = partialschur(shift_and_invert(A,-s*ε,σ),nev=2);

# decomp, history = partialschur(shift_and_invert(A,σ),nev=5);
k²,ψ = partialeigen(decomp,σ)
# display(sort(sqrt.((k)),by=real))
display(sqrt.(k²))
# p = reshape(ψ[:,1],N...)
# plot(abs2.(ψ[:,:]))
##
heatmap(abs2.(reshape(ψ[:,2],N...))')

##
using ScalarFDFD
##

∂Ω = [-1  -1
      +1  +1]
bnd = Boundary(∂Ω=∂Ω, bc=:d, bl=:pml, bl_depth=.3)    # boundary object
# bnd = Boundary(r=3, bl=bl, bl_depth=bl_depth)

dx = .01
dis = Discretization([dx,dx],Cartesian())
##
sys = System(Circle(π/6,0,0,Dict(:n1=>3)), Universe(1))
##
sim = Simulation(bnd=bnd, sys=sys, dis=dis)
plot(sim)

##

k, ψ = eig_kl(sim, 10, 1, 0, 0)
plot(sim, ψ, by=abs2, truncate=false)

##

plot(sim, ψ)
##

# sys = add_pc_waveguide(sys; width=1, direction=:vertical, x0=0, y0=0)
# sys = add_pc_waveguide(sys; width=1, direction=:horizontal, x0=0, y0=0)
#
# print(sys)
#
# sct = Scattering([Channels(1,1), Channels(2,1)])
#
# sim = Simulation(;sys=sys, bnd=bnd, dis=dis, sct=sct)
