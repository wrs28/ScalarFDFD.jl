# load ScalarFDFD on all process (@everywhere necessary b/c we will use parallel features)
@everywhere using ScalarFDFD

# define boundary of computation domain Ω:
∂Ω = [-8  -8
       9  9]

# and associate boundary conditions with each.
# Here all share the same, so we just specify a scalar.
bc = :d         # boundary condition (in this case Dirichlet)
bl = :pml       # boundary layer (in this case PML)
bl_depth = 3    # boundary layer depth

# now specify the lattice spacing (i.e. resolution)
dx = .05

# now we build the photonic crystal.
# first initialize the lattice and the pc object
pc = build_pc_domain(Bravais(a=1, b=1), Regions())
# then add two pillars, whose positions are referenced to the southwest corner of the cell
pc = add_circle_to_pc(pc, :sw; R=.1, n₁ = 3, x0=.3, y0=.3)
pc = add_circle_to_pc(pc, :sw; R=.20, n₁ = 3, x0=.7, y0=.7)

# now take all the definitions and build the various objects necessary for the simulation
bnd = Boundary(∂Ω=∂Ω, bc=bc, bl=bl, bl_depth=bl_depth)    # boundary object
dis = Discretization(dx)                                  # discretization object
sys = System(pc)                                          # system object

# now add waveguides. couldn't do this earlier because the different regions didn't
# know about each other until we built the "System" object.
# first we add the vertical waveguide, then the horizontal, each 1-unit cell wide
sys = add_pc_waveguide(sys; width=1, direction=:vertical)
sys = add_pc_waveguide(sys; width=1, direction=:horizontal)

# now print the results to see which waveguides correspond to which side
# this is in fact a deterministic recipe, but why bother? we see with the order
# in which we created these it goes 1: bottom, 2: top, 3: left, 4: right.
print(sys)

# specify what we define our channels to be. we could have done this before, but
# this way we know which waveguides we are referring to. btw the specification is
# Channel(waveguide_number, quantum_number). In this case, quantum number means
# transverse mode number. none of this is necessary for eigenproblems, only for
# scattering problems.
sct = Scattering([Channels(1,1), Channels(2,1)])

# now we are ready to put the whole simulation together into a neat package!
sim = Simulation(;sys=sys, bnd=bnd, dis=dis, sct=sct)

# to visualize the structure we just made,
plot(sim)
savefig("test1/pc_sim1.pdf")

# let's see what frequencies we should be operating at by looking at the band
# structure is. Note we specify parallel to be true, but in face it defaults to
# checking whether multiple processes are active and responding accordingly
bands, gaps = band_structure(sim; waveguide=1, zone=1, parallel=true)

# let's plot and save the result. all relevant plots are easiest if the first
# argument is the simulation. the order and number of the rest don't really matter
plot(sim, gaps, bands)
savefig("test1/pc_band_structure1.pdf")

# now build the dispersion curves for each waveguide. this is not necessary to
# do explicitly, as the first time it is needed for a scattering calculation
# it checks to see if it has already been done, and, if not, it does it and save
# the result to `sim`.
# but doing it this way allows much more control, though we are not using that here.
build_dispersions!(sim)

# to take a look at the dispersion curve, we just
plot(sim, sim.sct) #though plot(sim,sim) would have done just the same
savefig("test1/pc_wg_dispersion1.pdf")

# all this tells us that we can meaningfully expect to do a scattering calculation
# are ω=2.2 (say)
ψ_sct, ψ_tot = scattering(sim, 2.2, [1,0]) # the third argument is an array of amplitudes (one for each channel)

# let's look at the result!
plot(sim, ψ_tot)
savefig("test1/pc_scattering1.pdf")

# and make a cool animation
animate(wave(sim, ψ_tot, by=real), "test1/pc_scattering1.gif")


################################################################################

# let's do another example

∂Ω = [-7  -8*sqrt(3)/2
       14  8*sqrt(3)/2]

bc = :d
bl = :abs # this time we use the slightly more appropriate absorbing boundary layer for photonic crystals, but it must be longer
bl_depth = 4
bnd = Boundary(∂Ω=∂Ω, bc=bc, bl=bl, bl_depth=bl_depth)    # boundary object

dx = .05
dis = Discretization(dx)

# specify the angle of the second primitive vector. this choice corresponds to a
# hexagonal lattice
pc = build_pc_domain(Bravais(a=1, b=1, β=π/3), Regions())
# add a pillar to the middle of the cell
pc = add_circle_to_pc(pc, :center; R=.2, n₁ = 3)

# add a waveguide, this time it will terminate in the bulk
sys = System(pc)
sys = add_pc_waveguide(sys; width=1, direction=:west)

sct = Scattering([Channels(1,1)])

sim = Simulation(;sys=sys, bnd=bnd, dis=dis, sct=sct)

plot(sim)
savefig("test2/pc_sim2.pdf")

# this time let's map out the whole surface, but over the reduced zone
# (by making num_bloch an array instead of a number)
bands, gaps = band_structure(sim; waveguide=1, num_bloch=[21], zone=:reduced)
plot(sim, bands, gaps)
savefig("test2/pc_band_structure2.pdf")

bands, gaps = band_structure(sim; waveguide=1, num_bloch=51, zone=:reduced)

# same plotting as before now plots the surface
plot(sim, bands, gaps)
savefig("test2/pc_band_structure2_alt.pdf")

# we see there is a pretty wide gap between about 2 and 3. so let's venture a guess
# and scatter at ω = 2.5
# note this time we are not explicitly building the dispersion curves, as this
# is done during the scattering calculation. we can check after to see if our
# guess was a safe one
ψ_sct, ψ_tot = scattering(sim, 2.6, [1])

# let's again plot
plot(sim, ψ_tot, by=real)
savefig("test2/pc_scattering2.pdf")

# or perhaps you don't like that color scheme? or some other parameters?
plot(sim, ψ_tot, by=real, :solarized_light)
savefig("test2/pc_scattering2_alt.pdf")

# and let's make a movie, but with contours this time
animate(wave(sim, ψ_tot, :solarized_light; by=real, seriestype=:contour), "test2/pc_scattering2.gif")




################################################################################

# one last example

∂Ω = [-8  -8
       8   8]

bc = :d
bl = :abs
bl_depth = 4
bnd = Boundary(∂Ω=∂Ω, bc=bc, bl=bl, bl_depth=bl_depth)

dx = .05
dis = Discretization(dx)

pc = build_pc_domain(Bravais(a=1, b=1), Regions())
pc = add_circle_to_pc(pc, :center; R=.2, n₁ = 3)

sys = System(pc)
sys = add_pc_waveguide(sys; width=1, direction=:west)
sys = add_pc_waveguide(sys; width=1, direction=:south)

sct = Scattering([Channels(1,1)])

sim = Simulation(;sys=sys, bnd=bnd, dis=dis, sct=sct)

plot(sim)
savefig("test3/pc_sim3.pdf")

bands, gaps = band_structure(sim; waveguide=1)
plot(sim, bands, gaps)
savefig("test3/pc_band_structure3.pdf")

k, ψ = eig_k(sim, 2.5, :u, 1; direction=[1,-1])
println(k)
plot(sim,ψ,by=real)
savefig("test3/pc_uni_eig3.pdf")

ψ_sct, ψ_tot = scattering(sim, 2.5, [1])
plot(sim, ψ_tot, :default; by=real, seriestype=:contour)
savefig("test3/pc_scattering3.pdf")

# and let's make a movie, but with contours this time
animate(wave(sim, ψ_tot, :default; by=real), "test3/pc_scattering3.gif")
