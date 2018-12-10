################################################################################
### DOMAIN
################################################################################
"""
    domain = Domain(domain; :key1 => value1, :key2 => value2, ...)

New domain structure from old with modified fields
"""
function Domain(dom::Td;
    is_in_domain = is_in_domain,
    domain_params = domain_params,
    domain_type = domain_type,
    domain_ε = domain_ε,
    domain_F = domain_F,
    is_in_subdomain = is_in_subdomain,
    subdomain_params = subdomain_params,
    subdomain_type = subdomain_type,
    subdomain_ε = subdomain_ε,
    subdomain_F = subdomain_F,
    lattice = lattice,
    which_asymptote = which_asymptote,
    which_waveguide = which_waveguide
    ) where Td<:Domain
    return Domain(;is_in_domain = is_in_domain, domain_params = domain_params, domain_type = domain_type,
            domain_ε = domain_ε, domain_F = domain_F, is_in_subdomain = is_in_subdomain, subdomain_params = subdomain_params,
            subdomain_type = subdomain_type, subdomain_ε = subdomain_ε, subdomain_F = subdomain_F, lattice = lattice,
            which_asymptote = which_asymptote, which_waveguide = which_waveguide)
end


################################################################################
### SYSTEM
################################################################################
"""
    sys = System(domain1, domain2, ...)
"""
function System(args::Vararg{Td}) where Td<:Domain
    return System(vcat(args...))
end


"""
    sys = System(sys)

system object from system object
"""
function System(sys::Ts) where Ts<:System
    return System(sys.domains)
end


################################################################################
### DISCRETIZATION
################################################################################
"""
    dis = Discretization(dx, sub_pixel_num=DEFAULT_SUBSAMPLE_NUMBER, coordinate_system=:cart)

discretization object for Simulation

`dx` lattice spacing (scalar)

`sub_pixel_num` is the subsampling rate used in sub-pixel smoothing.
"""
function Discretization(dx::Real, coordinate_system::T=Cartesian(), sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER, origin=[0.,0.], N=[1,1], dN=[0 0;0 0]) where T<:CoordinateSystem
    return Discretization([dx,dx], coordinate_system, sub_pixel_num, origin, N, dN)
end
function Discretization(dx::Array, coordinate_system::T=Cartesian(), sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER) where T<:CoordinateSystem
    origin=[0.,0.]; N=[1,1]; dN=[0 0;0 0]
    return Discretization(dx, coordinate_system, sub_pixel_num, origin, N, dN)
end


"""
    dis = Discretization(dis; key1 => value1, key2 => value2...)

new discretization object from old, with modified fields.
"""
function Discretization(dis::Discretization; dx=dis.dx, coordinate_system=dis.coordinate_system, sub_pixel_num=dis.sub_pixel_num)
    return Discretization(dx, coordinate_system, sub_pixel_num, [0., 0.], [1,1], [0 0;0 0])
end



################################################################################
### BOUNDARY
################################################################################
"""
    bnd = Boundary(;∂Ω, bc=:d, bl=:none, bl_depth=0)

boundary object for Simulation.

`∂Ω` array of boundaries of computational domain. If not given, array of NaNs

`bc` scalar or array of boundary conditions

`bl` scalar or array of boundary layer types

`bl_depth` scalar or array of boundary layer depths
"""
function Boundary(;∂Ω=fill(NaN,2,2), bc=:d, bl=:none, bl_depth=0, r::Real=NaN)
    if !isnan(r)
        @assert typeof(bc)<:Symbol "type of bc is $(typeof(bc)). when radius is specified, bc must be Symbol"
        @assert typeof(bl)<:Symbol "type of bl is $(typeof(bl)). when radius is specified, bl must be Symbol"
        @assert typeof(bl_depth)<:Number "type of bl_depth is $(typeof(bl_depth)). when radius is specified, bl_depth must be Number"
        ∂Ω = [0   0
              r  2π]
        bc = [:d :p
              bc :p]
        bl = [:none :none
              bl    :none]
        bl_depth = [0        0
                    bl_depth 0]
    end
    return Boundary(∂Ω, bc, bl, bl_depth)
end
function Boundary(∂Ω,
        bc::Symbol, bl::Array{Symbol,2}, bl_depth::Array{T,2}) where T
    return Boundary(∂Ω, fill(bc,2,2), bl, bl_depth)
end
function Boundary(∂Ω, bc,
        bl::Symbol, bl_depth::Array{T,2}) where T
    return Boundary(∂Ω, bc, fill(bl,2,2), bl_depth)
end
function Boundary(∂Ω, bc, bl,
    bl_depth::Number)
    return Boundary(∂Ω, bc, bl, fill(bl_depth,2,2))
end


"""
    bnd = Boundary(bnd; :key1 => value1, :key2 => value2, ...)

new boundary object from old, with modified fields.
"""
function Boundary(bnd::Boundary; ∂Ω=bnd.∂Ω, bc=bnd.bc,
            bl=bnd.bl, bl_depth=bnd.bl_depth)
    return Boundary(∂Ω, bc, bl, bl_depth)
end



################################################################################
### CHANNELS
################################################################################
"""
    channel = Channels(chn; key1 => value1, key2 => value2, ...)

new channel object from old with modified fields.
"""
function Channels(chn::Channels; waveguide=chn.waveguide, quantum_number=chn.quantum_number)
    return Channels(waveguide, quantum_number)
end



################################################################################
### SCATTERING
################################################################################
function Scattering(channel::Channels)
    return Scattering([channel])
end
"""
    sct = Scattering(sct; :key1 => value1, :key2 => value2, ...)

new scattering object from old with modified fields.
"""
function Scattering(sct::Scattering; channels=sct.channels)
    return Scattering(channels)
end

"""
    sct = Scattering(channel1, channel2, ...)
"""
function Scattering(args::Vararg{Channels})
    return Scattering(vcat(args...))
end

################################################################################
### TWOLEVELSYSTEM
################################################################################
"""
    tls = TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...)

new tls object from old, with modified fields
"""
function TwoLevelSystem(tls::TwoLevelSystem; D₀=tls.D₀, k₀=tls.k₀, γp=tls.γp, ω₀=tls.ω₀)
        TwoLevelSystem(tls.D₀, tls.k₀, tls.γp, tls.ω₀)
end


################################################################################
### SIMULATION
################################################################################
"""
    sim = Simulation(domain; bnd=Boundary(bc=:p), dis, sct, tls, Na=1, Nb=1)

Simulation from a single domain.

Intended for case where `domain` has a finite lattice structure, in which case
    it defines the simulation to be `Na`×`Nb` unit cells.

Care must be taken with boundary layer specifications. It is recommended that
`bnd` be left unspecified.
"""
function Simulation(dom::Domain; bnd::Boundary=Boundary(bc=:p), dis::Discretization, sct::Scattering=Scattering(), tls::TwoLevelSystem=TwoLevelSystem(), Na=1, Nb=1)
    lattice = dom.lattice
    bnd = deepcopy(bnd)
    if !isinf(lattice.a)
        bnd.bc[:,1] .= :p
        bnd.∂Ω[1,1] = lattice.x0
        bnd.∂Ω[2,1] = lattice.x0+Na*lattice.a*(cos(lattice.α)-sin(lattice.α)*cot(lattice.β))
    end
    if !isinf(lattice.b)
        bnd.bc[:,2] .= :p
        bnd.∂Ω[1,2] = lattice.y0
        bnd.∂Ω[2,2] = lattice.y0+Nb*lattice.b*sin(lattice.β)
    end
    return deepcopy(Simulation(sys=System(dom); bnd=Boundary(bnd), dis=dis, tls=tls, lat=BravaisLattice(lattice; :a=>Na*lattice.a, :b=>Nb*lattice.b)))
end


"""
    sim = Simulation(domain, sim; Na=1, Nb=1)

Build simulation out of a single `domain`, with fields other than `bnd` the same
as in `sim`, except it is `Na`×`Nb` unit cells.

Intended for case where `domain` has a finite lattice structure.
"""
function Simulation(dom::Domain, sim::Simulation; Na=1, Nb=1)
    return deepcopy(Simulation(dom; dis=sim.dis, tls=sim.tls, Na=Na, Nb=Nb))
end


"""
    sim = Simulation(domain_num, sim; Na=1, Nb=1)

Build simulation out of domain number `domain_num`, with fields other than `bnd`
the same as in `sim`, except it is `Na`×`Nb` unit cells.

Intended for case where specified domain has a finite lattice structure.
"""
function Simulation(domain_num::Int, sim::Simulation; Na=1, Nb=1)
    return deepcopy(Simulation(sim.sys.domains[domain_num]; dis=sim.dis, tls=sim.tls, Na=Na, Nb=Nb))
end


"""
    sim = Simulation(sim; key1 => val1, key2 => val2, ...)

new simulation object from old, with modified fields.
"""
function Simulation(sim::Simulation; sys=sim.sys, bnd=sim.bnd, dis=sim.dis, sct=sim.sct, tls=sim.tls, lat=sim.lat, disp_opt=false)
    return deepcopy(Simulation(sys=System(sys), bnd=Boundary(bnd), dis=Discretization(dis),
            sct=Scattering(sct), tls=TwoLevelSystem(tls), lat=BravaisLattice(lat), disp_opt=disp_opt))
end


"""
    isCartesian(sim)
"""
function DifferentialOperators.isCartesian(sim::Simulation)
    return isCartesian(sim.dis)
end
function DifferentialOperators.isCartesian(sim::Discretization)
    return isCartesian(dis.coordinate_system)
end


"""
    isPolar(sim)
"""
function DifferentialOperators.isPolar(sim::Simulation)
    return isPolar(sim.dis)
end
function DifferentialOperators.isPolar(dis::Discretization)
    return isPolar(dis.coordinate_system)
end
