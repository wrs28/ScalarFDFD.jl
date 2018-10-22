using ScalarFDFD

∂Ω = [-4  -4
       4   4]

bc = :d
bl = :pml
bl_depth = 1

bnd = Boundary(∂Ω, bc, bl, bl_depth)

dx = .005
dis = Discretization(dx)

sys = System(Domain())

sct = Scattering([Channels(1,1), Channels(2,2), Channels(3,3)])

sim = Simulation(;sys=sys, bnd=bnd, dis=dis, sct=sct)
sim = add_vertical_planar_waveguides(sim,; width=.3, index=2, x0=0)
sim = add_horizontal_planar_waveguides(sim,; width=.3, index=2, y0=0)

build_dispersions!(sim)

_, ψ1, _, H = scattering(sim, 14, [1, 0 , 0])
_, ψ2 = scattering(sim, 14., [0, 1 , 0], H=H)
_, ψ3 = scattering(sim, 14., [0, 0 , 1], H=H)

plot(sim, ψ1)
savefig("test4/planar1.pdf")

plot(sim, ψ2)
savefig("test4/planar2.pdf")

plot(sim, ψ3)
savefig("test4/planar3.pdf")

_, ψtogether = scattering(sim, 14., [1, 0 , 3], H=H)
animate(wave(sim, ψtogether, by=real), "test4/planar4.gif")
