#TODO: transverse_mode_1: give transverse mode definite sign so it doesn't change from call to call
#TODO: incident_pc, pc_transverse_field: higher gap-order modes
#TODO: higher order gaps in build_channels!
#TODO: fix total kluge in waveguide_dispersion which eliminates bleeding of bands into gap

################################################################################
### SOURCE SYNTHESIS
################################################################################
"""
    j, φ₊₋, φ₊, φ₋ = synthesize_source(sim, k, a)

`a[n]` complex amplitude associated with channel `n`
`j` source for scattering calculation
`φ₊[:,m]` is incident wave on waveguide `m`
`φ₋[:,m]` is outgoing wave on waveguide `m`
`φ₊₋` is outgoing + incident wave at disordered sites
"""
function synthesize_source(sim::Simulation, k, a)

    N = prod(sim.dis.N)
    W = length(sim.sct.waveguides_used)

    φ₊₋ = zeros(ComplexF64,N,W)
    φ₊ = zeros(ComplexF64,N,W)  # incident
    φ₋ = zeros(ComplexF64,N,W)  # outgoing

    φt₊₋ = zeros(ComplexF64,N) # incident
    φt₊ = zeros(ComplexF64,N) # incident
    φt₋ = zeros(ComplexF64,N) # outgoing

    for m ∈ 1:length(sim.sct.channels)
        ind = findall(sim.sct.channels[m].waveguide .== sim.sct.waveguides_used)
        φt₊₋, φt₊, φt₋ = ScalarFDFD.incident_mode(sim, k, m)
        φ₊₋[:,ind] += a[m]*φt₊₋
        φ₊[:,ind] += a[m]*φt₊
        φ₋[:,ind] += a[m]*φt₋
    end

    k² = k^2
    k²δ = Array{ComplexF64}(undef,N)
    j = zeros(ComplexF64,N)
    for w ∈ 1:W
        φ₊₋[:,w] = (sim.sys.ε[:] .!== sim.sct.ε₀[w][:]).*φ₊₋[:,w]
        φ₊[:,w] = (sim.sys.ε[:] .== sim.sct.ε₀[w][:]).*φ₊[:,w]
        φ₋[:,w] = (sim.sys.ε[:] .== sim.sct.ε₀[w][:]).*φ₋[:,w]
        k²δ[:] = k²*(sim.sys.ε[:]-sim.sct.ε₀[w][:])
        j += -k²δ.*φ₊₋[:,w]
    end

    return j, sum(φ₊₋, dims=2), sum(φ₊, dims=2), sum(φ₋, dims=2)
end


"""
    φ₊₋, φ₊, φ₋ = incident_mode(sim, k, m)
"""
function incident_mode(sim::Simulation, k, m)
    if isempty(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, sim.sct.channels[m].waveguide)])
        φ₊₋, φ₊, φ₋ = incident_free_space_mode(sim, k, m)
    elseif ScalarFDFD.ispc(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, sim.sct.channels[m].waveguide)][end])
        φ₊₋, φ₊, φ₋ = incident_pc_mode(sim, k, m)
    elseif ishalfspace(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, sim.sct.channels[m].waveguide)][1])
        φ₊₋, φ₊, φ₋ = ScalarFDFD.incident_halfspace_mode(sim, k, m)
    elseif  ScalarFDFD.isplanar(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, sim.sct.channels[m].waveguide)][1])
        φ₊₋, φ₊, φ₋ = ScalarFDFD.incident_planar_mode(sim, k, m)
    else
        throw(ArgumentError("unrecognized waveguide type"))
    end

    return φ₊₋, φ₊, φ₋
end


################################################################################
### FREE SPACE MODES (ANGULAR MOMENTUM)
################################################################################
"""
    φ₊₋, φ₊, φ₋ = incident_free_space_mode(sim, k, m)
"""
function incident_free_space_mode(sim::Simulation, k, m)

    x = sim.dis.x[1]
    y = sim.dis.x[2]

    R = hypot.(x,y)
    θ = atan.(x,y)

    φ₊ = exp(1im*m*θ).*hankelh2(m,k*R)/2
    φ₋ = exp(1im*m*θ).*hankelh1(m,k*R)/2
    φ₊₋ = exp(1im*m*θ).*besselj(m,k*R)/4

    return φ₊₋, φ₊, φ₋
end

################################################################################
### PHOTONIC CRYSTAL WAVEGUIDE MODES
################################################################################
"""
    φ₊₋, φ₊, φ₋, φ_interpolation = incident_pc_mode(sim, k, m)
"""
function incident_pc_mode(sim::Simulation, k, m)

    prop_const, utp, wg_sim = ScalarFDFD.pc_transverse_field(sim, k, m)

    X = repeat(sim.dis.x[1],1,sim.dis.N[2])
    Y = repeat(sim.dis.x[2],sim.dis.N[1],1)
    φ = utp.(X,Y)

    if ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:left]
        v = wg_sim.lat.v1
        β = +prop_const[1]
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:right]
        v = wg_sim.lat.v1
        β = -prop_const[1]
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:bottom]
        v = wg_sim.lat.v2
        β = +prop_const[1]
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:top]
        v = wg_sim.lat.v2
        β = -prop_const[1]
    end

    φ₊ = φ.*exp.(+1im*β*(v[1]*X+v[2]*Y))/sqrt(abs(real(β)))
    φ₋ = zeros(ComplexF64, size(φ₊))
    φ₊₋ = φ₊ + φ₋

    return φ₊₋[:], φ₊[:], φ₋[:]
end


"""
    pc_transverse_field(sim, k, m)
"""
function pc_transverse_field(sim::Simulation, k, m)

    waveguide = sim.sct.channels[m].waveguide
    wg_sim = extract_waveguide_simulation(sim, waveguide)
    if sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, waveguide)][1].which_asymptote ∈ [:top, :bottom]
        β = (optimize(x->(sim.sct.channels[m].dispersion[1](x)-float(k))^2, 0, 2π/2wg_sim.lat.b).minimizer[1])::Float64
        k_new, ψ = eig_k(wg_sim, k, 1; kb=β)
        ψ = ψ.*exp.(-complex.(0,wg_sim.dis.x[2])*β .+ 0wg_sim.dis.x[1])[:]
    else
        β = (optimize(x->(sim.sct.channels[m].dispersion[1](x)-float(k))^2, 0, 2π/2wg_sim.lat.a).minimizer[1])::Float64
        k_new, ψ = eig_k(wg_sim, k, 1; ka=β)
        ψ = ψ.*exp.(-complex.(0,wg_sim.dis.x[1])*β .+ 0wg_sim.dis.x[2])[:]
    end
    if !isapprox(k_new[1],k; atol=1e-2)
        @warn "computed k $(k_new[1]) not consistent with precomputed dispersion $k."
    end
    ψ = ψ*exp(-complex(0,angle(ψ[findmax(abs.(ψ))[2]])))/sqrt(quadrature(wg_sim, abs2.(ψ[:])))
    xs = wg_sim.dis.x[1][1] .+ (1:wg_sim.dis.N[1])*wg_sim.dis.dx[1]
    ys = wg_sim.dis.x[2][1] .+ (1:wg_sim.dis.N[2])*wg_sim.dis.dx[2]
    utp = CubicSplineInterpolation((xs,ys),reshape(ψ[:,1],wg_sim.dis.N[1],wg_sim.dis.N[2]), bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    return β, utp, wg_sim
end


################################################################################
### HALFSPACE WAVEGUIDES MODES
################################################################################
"""
    φ₊₋, φ₊, φ₋ = incident_halfspace_mode(sim, k, m)
"""
function incident_halfspace_mode(sim::Simulation, k, m)

    c, prop_const = ScalarFDFD.halfspace_transverse_field(sim, k, m)

    if sim.dis.N[1]==1
        transverse_dim = 1
        longitudinal_dim = 2
    elseif sim.dis.N[2]==1
        transverse_dim = 2
        longitudinal_dim = 1
    else
        throw(ArgumentError("halfspace incident mode will not work for 2d"))
    end

    xs = sim.dis.x[longitudinal_dim][1] .+ (1:sim.dis.N[longitudinal_dim])*sim.dis.dx[longitudinal_dim]

    if ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:left]
        β = +prop_const
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:right]
        β = -prop_const
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:bottom]
        β = +prop_const
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:top]
        β = -prop_const
    end

    if issubset(sim.bnd.bl[:,longitudinal_dim],[:pml_out, :abs_out])
        φ₊ = exp.(1im*β*xs)/sqrt(abs(real(β)))
        φ₋ = zeros(ComplexF64, size(φ₊))
    elseif issubset([:none],sim.bnd.bl[:,longitudinal_dim])
        x0 = sim.bnd.bc[2,longitudinal_dim]
        if sim.bnd.bc[2,longitudinal_dim] == :d
            φ₊ = exp.(1im*β*(xs-x0))/sqrt(abs(real(β)))
            φ₋ = -1im*exp.(-1im*β*(xs-x0))/sqrt(abs(real(β)))/2
        elseif sim.bnd.bc[2,longitudinal_dim] == :n
            φ₊ = exp.(1im*β*(xs-x0))/sqrt(abs(real(β)))
            φ₋ = +1im*exp.(-1im*β*(xs-x0))/sqrt(abs(real(β)))/2
        else
            throw(ErrorException("for now one-sided scattering can only be done from the left or bottom"))
        end
    end
    φ₊₋ = φ₊ + φ₋

    return φ₊₋[:], φ₊[:], φ₋[:]
end


"""
    1, β = halfspace_transverse_field(sim, k, m)
"""
function halfspace_transverse_field(sim::Simulation, k, m)
    prop_const = (optimize(x->(sim.sct.channels[m].dispersion[1](x[1])-float(k))^2, [float(k)], BFGS()).minimizer[1])::Float64
    return 1, prop_const
end


################################################################################
### PLANAR WAVEGUIDE MODES
################################################################################
"""
    incident_planar_mode(sim, k, m)
"""
function incident_planar_mode(sim::Simulation, k, m)

    β, u_itp, longitudinal_dim, transverse_dim = ScalarFDFD.planar_transverse_field(sim, k, m)

    if k > β
        @warn "mode $(sim.sct.channels[m].quantum_number) on waveguide $(sim.sct.channels[m].waveguide) is not guided at frequency $k"
    end

    waveguide = sim.sct.channels[m].waveguide
    asymptote = ScalarFDFD.get_asymptote(sim, waveguide)

    if issubset(sim.bnd.bl[:,longitudinal_dim],[:pml_out, :abs_out])
        if asymptote ∈ [:left, :bottom]
            φ₊ = u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,+β)*sim.dis.x[longitudinal_dim])/sqrt(β)
        elseif asymptote ∈ [:right, :top]
            φ₊ = u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,-β)*sim.dis.x[longitudinal_dim])/sqrt(β)
        else
            throw(ErrorException("invalid aymptote $(asymptote)"))
        end
        φ₋ = zeros(ComplexF64, size(φ₊))
        φ₊₋ = copy(φ₊)
    elseif issubset([:none],sim.bnd.bl[:,longitudinal_dim])
        if asymptote ∈ [:left, :bottom]
            if sim.bnd.bc[2,longitudinal_dim] == :d
                φ₊ = +u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,+β)*(sim.dis.x[longitudinal_dim] .- sim.dis.x[longitudinal_dim][end]))/sqrt(β)
                φ₋ = -u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,-β)*(sim.dis.x[longitudinal_dim] .- sim.dis.x[longitudinal_dim][end]))/sqrt(β)
            elseif sim.bnd.bc[2,longitudinal_dim] == :n
                φ₊ = +u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,+β)*(sim.dis.x[longitudinal_dim] .- sim.dis.x[longitudinal_dim][end]))/sqrt(β)
                φ₋ = +u_itp.(sim.dis.x[transverse_dim]).*exp.(complex(0,-β)*(sim.dis.x[longitudinal_dim] .- sim.dis.x[longitudinal_dim][end]))/sqrt(β)
            end
        else
            throw(ErrorException("for now one-sided scattering can only be done from the left or bottom"))
        end
    else
        throw(ArgumentError("invalid boundary layer $(sim.bnd.bl)"))
    end
    φ₊₋ = φ₊ + φ₋

    return φ₊₋[:], φ₊[:], φ₋[:]
end

#@code_warntype checked
"""
    β, u_itp, longitudinal_dim, transverse_dim = planar_transverse_field(sim, k, m; η_init=-2, num_modes_multiplier=3)
"""
function planar_transverse_field(sim::Simulation, k, m)

    β = (optimize(x->(sim.sct.channels[m].dispersion[1](x[1])-float(k))^2, [float(k)], BFGS()).minimizer[1])::Float64

    wg_sim = extract_waveguide_simulation(sim, sim.sct.channels[m].waveguide)

    k_new, ψ = ScalarFDFD.planar_kl(wg_sim, β, k, 1)
    ψ = ψ*exp(-complex(0,angle(ψ[findmax(abs.(ψ))[2]])))/sqrt(quadrature(wg_sim, abs2.(ψ[:])))
    if !isapprox(k_new[1],k; atol=1e-2)
        @warn "computed k $(k_new[1]) not consistent with precomputed dispersion $k."
    end

    if wg_sim.dis.N[1]==1
        transverse_dim = 2
        longitudinal_dim=1
    else
        transverse_dim=1
        longitudinal_dim=2
    end
    xs = sim.dis.x[transverse_dim][1] .+ (1:sim.dis.N[transverse_dim])*sim.dis.dx[transverse_dim]
    utp = CubicSplineInterpolation(xs,real(ψ[:,1]), extrapolation_bc=Line())

    return β, utp, longitudinal_dim, transverse_dim
end

################################################################################
### AUXILLIARIES
################################################################################
"""
    waveguide_domains = get_waveguide_domains(sim::Simulation, waveguide)
"""
function get_waveguide_domains(sim::Simulation, waveguide)
    return findall([sim.sys.domains[i].which_waveguide for i ∈ eachindex(sim.sys.domains)] .== waveguide)
end


"""
    which_asymptote = get_asymptote(sim, waveguide, wg_inds)
"""
function get_asymptote(sim::Simulation, waveguide, wg_inds=get_waveguide_domains(sim, waveguide))
    return sim.sys.domains[wg_inds[1]].which_asymptote
end


"""
    wg_sim = extract_waveguide_simulation(sim, waveguide)
"""
function extract_waveguide_simulation(sim::Simulation, waveguide::Int)

    wg_inds = ScalarFDFD.get_waveguide_domains(sim, waveguide)

    if isempty(sim.sys.domains[wg_inds])
        throw(ErrorException("haven't done free space waveguide extraction yet"))
    elseif all(isplanar.(sim.sys.domains[wg_inds]))
        wg_sim = extract_planar_waveguide_simulation(sim, waveguide)
    elseif all(ispc.(sim.sys.domains[wg_inds]))
        wg_sim = extract_pc_waveguide_simulation(sim, waveguide)
    elseif all(ishalfspace.(sim.sys.domains[wg_inds]))
        throw(ErrorException("haven't done halfplane waveguide extraction yet"))
    else
        throw(ErrorException("invalid waveguide type"))
    end

    return deepcopy(wg_sim)
end


"""
    wg_sim = extract_planar_waveguide_simulation(sim, waveguide)
"""
function extract_planar_waveguide_simulation(sim::Simulation, waveguide::Int)

        wg_inds = ScalarFDFD.get_waveguide_domains(sim, waveguide)
        which_asymptote = ScalarFDFD.get_asymptote(sim, waveguide)

        domains = Array{Domain}(undef,length(wg_inds))
        for i ∈ eachindex(wg_inds)
              domains[i] = Domain(sim.sys.domains[wg_inds[i]]; :which_asymptote => :none, :which_waveguide => 0)
        end
        waveguide = System(domains)

        sim = deepcopy(sim)

        if which_asymptote == :right
            sim.bnd.∂Ω[1,1] = sim.bnd.∂Ω_tr[2,1]
            sim.dis.dx[1] = Inf
        elseif which_asymptote == :left
            sim.bnd.∂Ω[1,1] = sim.bnd.∂Ω_tr[1,1]
            sim.dis.dx[1] = Inf
        elseif which_asymptote == :bottom
            sim.bnd.∂Ω[1,2] = sim.bnd.∂Ω_tr[1,2]
            sim.dis.dx[2] = Inf
        elseif which_asymptote == :top
            sim.bnd.∂Ω[1,2] = sim.bnd.∂Ω_tr[2,2]
            sim.dis.dx[2] = Inf
        else
            throw(ArgumentError("Invalid asymptotic region $(which_asymptote)."))
        end
        sim.bnd.bl[:] .= :none
        sim.bnd.bl_depth[:] .= 0
        sim.bnd.bc[:] .= :p

        wg_sim = Simulation(sys=waveguide, dis=sim.dis, bnd=Boundary(sim.bnd))

        return wg_sim
end


"""
    wg_sim = extract_pc_waveguide_simulation(sim, waveguide)
"""
function extract_pc_waveguide_simulation(sim::Simulation, waveguide::Int)

    sim = deepcopy(sim)

    bnd = sim.bnd
    sys = sim.sys

    wg_inds = ScalarFDFD.get_waveguide_domains(sim, waveguide)
    which_asymptote = ScalarFDFD.get_asymptote(sim, waveguide, wg_inds)

    domains = Array{Domain}(undef,length(wg_inds))
    for i ∈ eachindex(wg_inds)
          domains[i] = Domain(sys.domains[wg_inds[i]]; :which_asymptote => :none, :which_waveguide => 0)
    end
    waveguide = System(domains)

    wg_inds = findall([waveguide.domains[i].domain_type for i ∈ eachindex(waveguide.domains)] .== :pc_waveguide)
    pc_inds = findall([waveguide.domains[i].domain_type for i ∈ eachindex(waveguide.domains)] .== :pc_waveguide_background)
    if isempty(wg_inds)
        ind = 1
    else
        ind = wg_inds[1]
    end
    if which_asymptote == :right
        waveguide_lattice = Bravais(waveguide.domains[ind].lattice; :x0 => bnd.∂Ω_tr[2,1], :y0 => bnd.∂Ω[1,2], :b => (bnd.∂Ω[2,2]-bnd.∂Ω[1,2])/sin(waveguide.domains[ind].lattice.β))
        bnd.∂Ω[1,1] = bnd.∂Ω_tr[2,1]
        bnd.∂Ω[2,1] = bnd.∂Ω[1,1] + waveguide_lattice.a*waveguide_lattice.v1[1]
    elseif which_asymptote == :left
        waveguide_lattice = Bravais(waveguide.domains[ind].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[1,2], :b => (bnd.∂Ω[2,2]-bnd.∂Ω[1,2])/sin(waveguide.domains[ind].lattice.β))
        bnd.∂Ω[2,1] = bnd.∂Ω_tr[1,1]
        bnd.∂Ω[1,1] = bnd.∂Ω[2,1] - waveguide_lattice.a*waveguide_lattice.v1[1]
    elseif which_asymptote == :bottom
        waveguide_lattice = Bravais(waveguide.domains[ind].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[1,2], :a => (bnd.∂Ω[2,1]-bnd.∂Ω[1,1])/cos(waveguide.domains[ind].lattice.α))
        bnd.∂Ω[2,2] = bnd.∂Ω_tr[1,2]
        bnd.∂Ω[1,2] = bnd.∂Ω[2,2] - waveguide_lattice.b*waveguide_lattice.v2[2]
    elseif which_asymptote == :top
        waveguide_lattice = Bravais(waveguide.domains[ind].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[2,2], :a => (bnd.∂Ω[2,1]-bnd.∂Ω[1,1])/cos(waveguide.domains[ind].lattice.α))
        bnd.∂Ω[1,2] = bnd.∂Ω_tr[2,2]
        bnd.∂Ω[2,2] = bnd.∂Ω[1,2] + waveguide_lattice.b*waveguide_lattice.v2[2]
    end
    if which_asymptote ∈ [:left, :right]
        bnd.bc[:,1] .= :p
        bnd.bl[:,1] .= :none
        bnd.bl_depth[:,1] .= 0
    elseif which_asymptote ∈ [:bottom, :top]
        bnd.bc[:,2] .= :p
        bnd.bl[:,2] .= :none
        bnd.bl_depth[:,2] .= 0
    else
        throw(ArgumentError("Invalid asymptotic region $(which_asymptote)."))
    end

    bnd.bc[:] .= :p
    bnd.bl[:] .= :none
    bnd.bl_depth[:] .= 0
    wg_sim = Simulation(sys=waveguide, dis=sim.dis, bnd=Boundary(bnd), lat=waveguide_lattice)

    return wg_sim
end
