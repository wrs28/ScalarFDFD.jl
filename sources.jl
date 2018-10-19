#TODO: transverse_mode_1: give transverse mode definite sign so it doesn't change from call to call
#TODO: incident_pc, pc_transverse_field: higher gap-order modes
#TODO: extract_waveguide_sim, build_channels! for fibers and in 1-dim
#TODO: higher order gaps in build_channels!
#TODO: fix total kluge in waveguide_dispersion which eliminates bleeding of bands into gap
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
        φt₊₋, φt₊, φt₋ = incident_mode(sim, k, m)
        φ₊₋[:,sim.sct.channels[m].waveguide] += a[m]*φt₊₋
        φ₊[:,sim.sct.channels[m].waveguide] += a[m]*φt₊
        φ₋[:,sim.sct.channels[m].waveguide] += a[m]*φt₋
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
    if ispc(sim.sys.domains[get_waveguide_domains(sim, sim.sct.channels[m].waveguide)][end])
        return incident_mode_pc(sim, k, m)
    else
        return incident_mode_planar(sim, k, m)
    end
end

################################################################################
### PHOTONIC CRYSTAL WAVEGUIDES
################################################################################
"""
    φ₊₋, φ₊, φ₋, φ_interpolation = incident_mode_pc(sim, k, m)
"""
function incident_mode_pc(sim::Simulation, k, m)
    prop_const, ψ, wg_sim = pc_transverse_field(sim, k, m)

    xs = wg_sim.dis.x[1][1] .+ (1:wg_sim.dis.N[1])*wg_sim.dis.dx[1]#wg_sim.dis.x[1][1]:(wg_sim.dis.x[1][end]-wg_sim.dis.x[1][1])/wg_sim.dis.N[1]:wg_sim.dis.x[1][end]
    ys = wg_sim.dis.x[2][1] .+ (1:wg_sim.dis.N[2])*wg_sim.dis.dx[2]#wg_sim.dis.x[2][1]:(wg_sim.dis.x[2][end]-wg_sim.dis.x[2][1])/wg_sim.dis.N[2]:wg_sim.dis.x[2][end]
    utp = CubicSplineInterpolation((xs,ys),reshape(ψ[:,1],wg_sim.dis.N[1],wg_sim.dis.N[2]), bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    XY = bravais_coordinates_unit_cell.(sim.dis.x[1],sim.dis.x[2], Ref(wg_sim.lat))
    X = Array{Float64}(undef,size(XY))
    Y = Array{Float64}(undef,size(XY))
    for i ∈ eachindex(XY)
          X[i]=XY[i][1]
          Y[i]=XY[i][2]
    end
    φ = utp.(X,Y)

    XY = bravais_coordinates.(sim.dis.x[1],sim.dis.x[2], Ref(wg_sim.lat))
    X = Array{Float64}(undef,size(XY))
    Y = Array{Float64}(undef,size(XY))
    for i ∈ eachindex(XY)
          X[i]=XY[i][1]
          Y[i]=XY[i][2]
    end

    if get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:left]
        v = wg_sim.lat.v1
        β = +prop_const[1]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:right]
        v = wg_sim.lat.v1
        β = -prop_const[1]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:bottom]
        v = wg_sim.lat.v2
        β = +prop_const[1]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) ∈ [:top]
        v = wg_sim.lat.v2
        β = -prop_const[1]
    end

    φ₊ = φ.*exp.(1im*β*(v[1]*X+v[2]*Y))/sqrt(abs(real(β)))
    φ₋ = zeros(ComplexF64, size(φ₊))
    φ₊₋ = copy(φ₊)
    return φ₊₋[:], φ₊[:], φ₋[:], utp, prop_const
end


"""
    pc_transverse_field(sim, k, m)
"""
function pc_transverse_field(sim::Simulation, k, m)
    waveguide = sim.sct.channels[m].waveguide
    wg_sim = extract_waveguide_simulation(sim, waveguide)
    if sim.sys.domains[get_waveguide_domains(sim, waveguide)][1].which_asymptote ∈ [:top, :bottom]
        β = (optimize(x->(sim.sct.channels[m].dispersion[1](x)-k)^2, 0, 2π/2wg_sim.lat.b).minimizer[1])::Float64
        k_new, ψ = eig_k(wg_sim, k, 1; kb=β)
        ψ = ψ.*exp.(complex.(0,wg_sim.dis.x[2])*β .+ 0wg_sim.dis.x[1])[:]
    else
        β = (optimize(x->(sim.sct.channels[m].dispersion[1](x)-k)^2, 0, 2π/2wg_sim.lat.a).minimizer[1])::Float64
        k_new, ψ = eig_k(wg_sim, k, 1; ka=β)
        ψ = ψ.*exp.(complex.(0,wg_sim.dis.x[1])*β .+ 0wg_sim.dis.x[2])[:]
    end
    if !isapprox(k_new[1],k; atol=1e-2)
        @warn "computed k $(k_new[1]) not consistent with precomputed dispersion $k."
    end
    ψ = ψ*exp(-complex(0,angle(ψ[findmax(abs.(ψ))[2]])))/sqrt(quadrature(wg_sim, abs2.(ψ[:])))
    return β, ψ, wg_sim
end

################################################################################
### PLANAR WAVEGUIDES
################################################################################

function incident_mode_planar(sim::Simulation, k, m; η_init=-2)

    β, u_itp, collapsed_dim, full_dim = planar_transverse_field(sim::Simulation, k, m; η_init=η_init)
    side = sim.bnd.waveguides[sim.sct.channels[m].waveguide].side

    if issubset(sim.bnd.bl[:,collapsed_dim],[:pml_out, :abs_out])
        if side ∈ [:l, :b]
            φ₊ = u_itp(sim.dis.x[full_dim]).*exp.(complex(0,+β)*sim.dis.x[collapsed_dim])/sqrt(β)
        elseif side ∈ [:r, :t]
            φ₊ = u_itp(sim.dis.x[full_dim]).*exp.(complex(0,-β)*sim.dis.x[collapsed_dim])/sqrt(β)
        else
            throw(ErrorException("Invalid side $(side)"))
        end
        φ₋ = zeros(ComplexF64, size(φ₊))
        φ₊₋ = copy(φ₊)
    elseif issubset([:none],sim.bnd.bl[:,collapsed_dim])
        if side ∈ [:l, :b]
            if sim.bnd.bc[2,collapsed_dim] == :d
                φ₊ = +u_itp(sim.dis.x[full_dim]).*exp.(complex(0,+β)*(sim.dis.x[collapsed_dim] .- sim.dis.x[collapsed_dim][end]))/sqrt(β)
                φ₋ = -u_itp(sim.dis.x[full_dim]).*exp.(complex(0,-β)*(sim.dis.x[collapsed_dim] .- sim.dis.x[collapsed_dim][end]))/sqrt(β)
            elseif sim.bnd.bc[2,collapsed_dim] == :n
                φ₊ = +u_itp(sim.dis.x[full_dim]).*exp.(complex(0,+β)*(sim.dis.x[collapsed_dim] .- sim.dis.x[collapsed_dim][end]))/sqrt(β)
                φ₋ = +u_itp(sim.dis.x[full_dim]).*exp.(complex(0,-β)*(sim.dis.x[collapsed_dim] .- sim.dis.x[collapsed_dim][end]))/sqrt(β)
            else
                throw(ErrorException("Haven't done half-open yet"))
            end
        end
        φ₊₋ = φ₊ + φ₋
    end
    return φ₊₋[:], φ₊[:], φ₋[:]
end

function planar_transverse_field(sim::Simulation, k, m; η_init=-2)

    sim = deepcopy(sim)

    side = sim.bnd.waveguides[sim.sct.channels[m].waveguide].side
    waveguide = sim.sct.channels[m].waveguide
    quantum_number = sim.sct.channels[m].quantum_number
    domain = sim.bnd.waveguide_domains[sim.sct.channels[m].waveguide]
    for i ∈ eachindex(domain.F_val)
        domain.F_val[i] = 1
    end
    sys = System(domain)

    if side == :l
        collapsed_dim = 1
        full_dim = 2
        ∂Ω = sim.bnd.∂Ω[1,collapsed_dim]
    elseif side == :r
        collapsed_dim = 1
        full_dim = 2
        ∂Ω = sim.bnd.∂Ω[2,collapsed_dim]
    elseif side == :b
        collapsed_dim = 2
        full_dim = 1
        ∂Ω = sim.bnd.∂Ω[1,collapsed_dim]
    elseif side == :t
        collapsed_dim = 2
        full_dim = 1
        ∂Ω = sim.bnd.∂Ω[2,collapsed_dim]
    else
        throw(ArgumentError("Invalid side specification $(side) for channel $(m)"))
    end

    sim.bnd.∂Ω_tr[1,collapsed_dim] = ∂Ω-.01
    sim.bnd.∂Ω_tr[2,collapsed_dim] = ∂Ω+.01
    sim.bnd.bl[:,collapsed_dim] .= :none
    sim.bnd.bl_depth[:,collapsed_dim] .= 0
    bnd = Boundary(sim.bnd; :waveguides => Waveguide[])

    sim.dis.N_tr[collapsed_dim] = 1
    sim.dis.N_tr[full_dim] = TRANSVERSE_MODE_RESOLUTION_MULTIPLIER::Int*sim.dis.N_tr[full_dim]

    sct = Scattering()

    sim = Simulation(sim; :bnd => bnd, :sys => sys, :sct => sct, :tls => TwoLevelSystem())

    β, u = eig_βl(sim, k, quantum_number; η_init=η_init)

    xs = sim.dis.x[full_dim][1]:sim.dis.dx[full_dim]:sim.dis.x[full_dim][end]
    u_itp = CubicSplineInterpolation(xs, u[:,quantum_number])

    return β[quantum_number], u_itp, collapsed_dim, full_dim
end

################################################################################
### DISPERSION BUILDING
################################################################################
"""
    build_dispersion!(sim; num_wg_bands=15, num_bloch=17, num_free_bands=2)
"""
function build_dispersion!(sim::Simulation;
    num_wg_bands_multiplier=1.8,
    num_bloch=17, num_free_bands=2, parallel=false)

    num_wg = length(sim.sys.waveguides)
    for i ∈ eachindex(sim.sys.waveguides)
        ch_inds = findall([sim.sct.channels[j].waveguide==sim.sys.waveguides[i] for j ∈ eachindex(sim.sct.channels)])
        if any(isempty.([sim.sct.channels[j].dispersion for j ∈ ch_inds]))
            println("Building dispersion for waveguide $i")
            wg_bands, bands, gaps, ks = waveguide_dispersion(sim, sim.sys.waveguides[i]; num_wg_bands_multiplier=num_wg_bands_multiplier, num_bloch=num_bloch, num_free_bands=num_free_bands, parallel=parallel)
            if isempty(wg_bands)
                throw(ErrorException("no guided modes found. Either there are none in the first bandgap, or increase num_wg_bands from $num_wg_bands."))
            end
            for j ∈ ch_inds
                if isempty(sim.sct.channels[j].dispersion)
                    push!(sim.sct.channels[j].dispersion, wg_bands[sim.sct.channels[j].quantum_number])
                    push!(sim.sct.channels[j].gaps, [gaps[1][1], gaps[2][1]])
                end
            end
        end
    end
    return nothing
end


"""
    remove_dispersion!(sim)

removes dispersion curves from simulation
"""
function remove_dispersion!(sim::Simulation)
    for i ∈ 1:length(sim.sct.channels)
        for j ∈ 1:length(sim.sct.channels[i].dispersion)
            pop!(sim.sct.channels[i].dispersion)
        end
    end
    return nothing
end


"""
    wg_bands, gaps = waveguide_dispersion(sim, waveguide; num_free_bands=2, num_wg_bands=17, num_bloch=17, interpolation=:cubic)
"""
function waveguide_dispersion(sim::Simulation, waveguide::Int;
    num_wg_bands_multiplier=1.8,
    num_free_bands=2, num_bloch=17, parallel=false, interpolation=:cubic)

    _, gaps = band_structure(sim; waveguide=waveguide, num_bands=num_free_bands, num_bloch=num_bloch, interpolation=interpolation, parallel=parallel)
    wg_dispersion, (bands, ks) = waveguide_dispersion(sim, waveguide, gaps; num_wg_bands_multiplier=num_wg_bands_multiplier, num_bloch=num_bloch, interpolation=interpolation, num_free_bands=num_free_bands, parallel=parallel)
    return wg_dispersion, bands, gaps, ks
end


"""
    wg_bands = waveguide_dispersion(sim::Simulation, waveguide, gaps; num_free_bands=2, num_wg_bands=15, num_bloch=17, interpolation=:cubic)
"""
function waveguide_dispersion(sim::Simulation, waveguide::Int, gaps::Tuple;
    num_wg_bands_multiplier=1.8,
    num_free_bands=2, num_bloch=17, parallel=false, interpolation=:cubic)

    wg_sim = extract_waveguide_simulation(sim, waveguide)
    num_wg_bands=round(Int,num_wg_bands_multiplier*max(wg_sim.lat.a/wg_sim.lat.b,wg_sim.lat.b/wg_sim.lat.a))
    pc_sim = Simulation(wg_sim.sys.domains[2]; bnd=Boundary(bc=:p), dis=wg_sim.dis)
    if wg_sim.lat.a < wg_sim.lat.b
        ka_bloch = 0:2π/wg_sim.lat.a/(num_bloch-1):2π/wg_sim.lat.a
        kb_bloch = fill(0,length(ka_bloch))
        k_bloch = ka_bloch
        ks = ka_bloch[1]:(ka_bloch[2]-ka_bloch[1]):ka_bloch[end]
    else
        kb_bloch = 0:2π/wg_sim.lat.b/(num_bloch-1):2π/wg_sim.lat.b
        ka_bloch = fill(0,length(kb_bloch))
        k_bloch = kb_bloch
        ks = kb_bloch[1]:(kb_bloch[2]-kb_bloch[1]):kb_bloch[end]
    end

    pg = Progress(Int(1e5), 1e5)
    free_bands, _, _ = band_structure(pc_sim, (ka_bloch, kb_bloch); num_bands=num_free_bands, parallel=parallel, interpolation=interpolation, pg=pg)
    _, _, bands = band_structure(wg_sim, (ka_bloch, kb_bloch); num_bands=num_wg_bands, parallel=parallel, calc_type="waveguide dispersion ")

    mode_bool = fill(false, num_wg_bands)
    for i ∈ eachindex(bands)
        mode_bool[i] = any(gaps[1][1]+.005 .< bands[i].(1:.01:num_bloch) .< gaps[2][1]-.005)
    end
    mode_inds = findall(mode_bool)

    wg_bands = Array{AbstractInterpolation}(undef,length(mode_inds))
    ks = k_bloch[1]:(k_bloch[2]-k_bloch[1]):k_bloch[end]
    for j in 1:length(mode_inds)
        wg_bands[j] = scale(bands[mode_inds[j]], ks)
    end
    return wg_bands, (free_bands, ks)
end

################################################################################
### AUXILLIARIES
################################################################################
"""
    ispc(dom::Domain)
"""
function ispc(dom::Domain)
    return dom.domain_type ∈ [:pc, :pc_waveguide]
end


"""
    wg_sim = extract_waveguide_simulation(sim, waveguide)
"""
function extract_waveguide_simulation(sim::Simulation, waveguide)

    old_system_size = sim.bnd.∂Ω_tr[2,:]-sim.bnd.∂Ω_tr[1,:]
    sim = deepcopy(sim)

    bnd = sim.bnd
    sys = sim.sys

    wg_inds = get_waveguide_domains(sim, waveguide)
    which_asymptote = get_asymptote(sim, waveguide, wg_inds)

    domains = Array{Domain}(undef,length(wg_inds))
    for i ∈ eachindex(wg_inds)
          domains[i] = Domain(sys.domains[wg_inds[i]]; :which_asymptote => :none, :which_waveguide => 0)
    end
    waveguide = System(domains)

    pc_inds = findall([waveguide.domains[i].domain_type for i ∈ eachindex(waveguide.domains)] .== :pc_waveguide)
    defect_inds = findall([waveguide.domains[i].domain_type for i ∈ eachindex(waveguide.domains)] .== :pc_waveguide_defect)
    if which_asymptote == :right
        waveguide_lattice = Bravais(waveguide.domains[pc_inds[1]].lattice; :x0 => bnd.∂Ω_tr[2,1], :y0 => bnd.∂Ω[1,2], :b => bnd.∂Ω[2,2]-bnd.∂Ω[1,2])
        bnd.∂Ω[1,1] = bnd.∂Ω_tr[2,1]
        bnd.∂Ω[2,1] = bnd.∂Ω[1,1] + waveguide_lattice.a*waveguide_lattice.v1[1]
    elseif which_asymptote == :left
        waveguide_lattice = Bravais(waveguide.domains[pc_inds[1]].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[1,2], :b => bnd.∂Ω[2,2]-bnd.∂Ω[1,2])
        bnd.∂Ω[2,1] = bnd.∂Ω_tr[1,1]
        bnd.∂Ω[1,1] = bnd.∂Ω_tr[2,1] - waveguide_lattice.a*waveguide_lattice.v1[1]
    elseif which_asymptote == :bottom
        waveguide_lattice = Bravais(waveguide.domains[pc_inds[1]].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[1,2], :a => bnd.∂Ω[2,1]-bnd.∂Ω[1,1])
        bnd.∂Ω[2,2] = bnd.∂Ω_tr[1,2]
        bnd.∂Ω[1,2] = bnd.∂Ω_tr[2,2] - waveguide_lattice.b*waveguide_lattice.v2[2]
    elseif which_asymptote == :top
        waveguide_lattice = Bravais(waveguide.domains[pc_inds[1]].lattice; :x0 => bnd.∂Ω[1,1], :y0 => bnd.∂Ω[2,2], :a => bnd.∂Ω[2,1]-bnd.∂Ω[1,1])
        bnd.∂Ω[1,2] = bnd.∂Ω_tr[2,2]
        bnd.∂Ω[2,2] = bnd.∂Ω_tr[2,2] + waveguide_lattice.b*waveguide_lattice.v2[2]
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


"""
    get_waveguide_domains(sim::Simulation, waveguide)
"""
function get_waveguide_domains(sim::Simulation, waveguide)
    wg_inds = findall([sim.sys.domains[i].which_waveguide for i ∈ eachindex(sim.sys.domains)] .== waveguide)
    return wg_inds
end

"""
    get_asymptote(sim, waveguide, wg_inds)
"""
function get_asymptote(sim::Simulation, waveguide, wg_inds=get_waveguide_domains(sim, waveguide))
    which_asymptote = sim.sys.domains[wg_inds[1]].which_asymptote
    return which_asymptote
end



#
# """
#     eig_βl(sim, k, nβ; η_init)
#
# propagation constant and transverse field for planar waveguides
# """
# function eig_βl(sim::Simulation, k::Number, nβ::Int; η_init=0)
#
#     η, v = eig_cf(sim, k, 3nβ; η_init=η_init)
#
#     fulldim = findall(sim.dis.N_tr .> 1)[1]
#
#     if issubset(sim.bnd.bl[:,full_dim],[:pml_out, :abs_out])
#         η = η[real(η) .< -1]
#     elseif issubset(sim.bnd.bl[:,full_dim], [:none])
#         η = η[real(η) .< 0]
#     else
#         @warn "haven't considered a half-open case"
#     end
#     γ = k*sqrt.(η)
#     α = real(γ)
#     perm = findall(α .< PROPAGATION_CONSTANT_IMAGINARY_CUTOFF::Float64)
#     β = abs.(imag(γ[perm]))
#     if length(β) < nβ
#         throw(ErrorException("Requested quantum number $(nβ) not supported at frequency $(k) in waveguide $(waveguide)"))
#     end
#     perms = sortperm(β, rev=true)
#
#     u = Array{Float64}(undef, sim.dis.N[full_dim], length(β))
#     for i ∈ 1:length(β)
#         u[:,i] = flipsign.(real(v[:,perm[perms[i]]]), real(v[findmax(abs.(real(v[:,perm[perms[i]]])))[2]]) )
#     end
#
#     return β[perms], u
# end


################################################################################
#### PROPAGATION CONSTANT SOLVERS
################################################################################
"""
    β, u = eig_β(sim, k, nβ; β=k, radii=(β,.001), η_init=0, rank_tol=1e-8, Nq=50)
"""
function eig_β(sim::Simulation, k, nβ; β=k, radii=(β,.001), η_init=0, rank_tol=1e-8, Nq=50)

    collapsed_dim = findall(sim.dis.N_tr .== 1)
    if !isempty(collapsed_dim)
        β, u = eig_β(sim, k, nβ; η_init=η_init)
    else
        β, u = eig_βnl(sim, k, β, 3nβ, radii; Nq=Nq, rank_tol=rank_tol)
    end
end


"""
    eig_β(sim, k; β=k, η_init=0, β_avoid=[0], disp_opt=false, tol=.5, max_count=15, max_iter=50)
"""
function eig_β(sim::Simulation, k; β=k, η_init=0, β_avoid=[0],
    disp_opt=false, tol=.5, max_count=15, max_iter=50)

    collapsed_dim = findall(sim.dis.N_tr .== 1)
    if !isempty(collapsed_dim)
        β, u = eig_β(sim, k, 1; η_init=η_init)
    else
        β, u = eig_βnl(sim, k, β; η_init=η_init, β_avoid=β_avoid,
            disp_opt=disp_opt, tol=tol, max_count=max_count, max_iter=max_iter)
    end
end



################################################################################
#### PROPAGATION CONSTANT SOLVERS
################################################################################
"""
    k, ψ, η, conv = eig_βnl(sim, k, β_init=k; F=[1], η_init=0, u_init=[], k_avoid=[0],
        disp_opt=false, tol=.5, max_count=15, max_iter=50)

eigenfrequency `k` with dispersion via root-finding of `η(k) - D₀γ(k)`.

`max_count`: the max number of CF states to compute at a time.
`max_iter`: maximum number of iterations to include in nonlinear solve
"""
function eig_βnl(sim::Simulation, k::Number, β_init::Number=k;
    F=[1], η_init=0, u_init=[], β_avoid=[0],
    disp_opt=false, tol=.5, max_count=15, max_iter=50)

    periodic_dim = findall(.!isinf.([sim.lat.a,sim.lat.b]))[1]

    if periodic_dim == 1
        symbol = :ka
    elseif periodic_dim == 2
        symbol = :kb
    else
        throw(ArgumentError("simulation is not a pc waveguide simulation"))
    end
    ηt, ut = eig_cf(sim, k, 1; u_init=u_init, symbol => β_init)

    γp = sim.tls.γp; k₀ = sim.tls.k₀; D₀ = sim.tls.D₀; dx = sim.dis.dx
    function f!(fvec,z)

        β = complex(z[1],z[2])

        flag = true
        count = 1
        M = 1
        ind = Int

        while flag
            η_temp, u_temp = eig_cf(sim, k, M; F=F, η_init=ηt[1], u_init=ut[:,1], symbol => β)
            overlap = zeros(Float64,M)
            for i ∈ eachindex(overlap)
                overlap[i] = abs(quadrature(sim,ut[:,1].*sim.sys.F[:].*F.*u_temp[:,i]))
            end
            max_overlap, ind = findmax(overlap)
            if max_overlap > (1-tol)
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
            elseif count < max_count && max_overlap ≤ (1-tol)
                M += 1
            else
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
                @warn "overlap less than tolerance, skipped to neighboring TCF."
            end
            count += 1
        end
        fvec[1] = real((ηt[1]-D₀*γ(sim,k))/prod(β .- β_avoid))
        fvec[2] = imag((ηt[1]-D₀*γ(sim,k))/prod(β .- β_avoid))
    end

    z = nlsolve(f!, [real(float(β_init)),imag(float(β_init))]; iterations=max_iter, show_trace=disp_opt)
    k = complex(z.zero[1], z.zero[2])
    conv = converged(z)

    ηt[1] = D₀*γ(sim,k)

    η, ψ = eig_cf(sim, k, 1; F=F, η_init=ηt[1], u_init = ut[:,1], symbol=>β)

    if !conv
        warn("no convergence for frequency $(k[1])")
    end
    return [β]::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}, η::Array{ComplexF64,1}, conv::Bool
end


"""
    β = eig_βnl(sim, k, nk, radii; F=[1], r_min=.01, Nq=100, rank_tol=1e-8, ka=0, kb=0)

Eigenfrequency with dispersion, using contour integration. BC's and BL's  set
by `sim.bnd.bc`, `sim.bnd.bl`

Contour is centered on `k`, `radii` = [real radius, imag radius].

`nk` is an upper bound on the number of eigenfrequencies contained in the contour.

`Nq` is the number of contour quadrature points.
"""
function eig_βnl(sim::Simulation, k::Number, β::Number, nβ::Int, radii; F=[1], r_min=.01, Nq=100, rank_tol=1e-8)

    ∇², ∇₁, ∇₂ = laplacian_components(sim, k)

    periodic_dim = findall(.!isinf.([sim.lat.a,sim.lat.b]))[1]
    if periodic_dim == 1
        ∇ = ∇₁
    elseif periodic_dim == 2
        ∇ = ∇₂
    else
        throw(ArgumentError("simulation is not a pc waveguide simulation"))
    end

    N = prod(sim.dis.N)
    A  = zeros(ComplexF64, N, nβ)
    A₀ = zeros(ComplexF64, N, nβ)
    A₁ = zeros(ComplexF64, N, nβ)
    M = rand(N, nβ)

    ϕ = 2π*range(0,length=Nq,step=1/Nq)
    B = β .+ complex.(radii[1]*cos.(ϕ), radii[2]*sin.(ϕ))

    ε=sim.sys.ε; F0=sim.sys.F; D₀=sim.tls.D₀; dx = sim.dis.dx; k² = k^2
    for i in 1:Nq

        β′ = B[i]
        β′² = sparse(1:N, 1:N, β′^2, N, N)

        if i > 1 && i < Nq
            dβ′ = (B[i+1] - B[i-1]  )/2
        elseif i == Nq
            dβ′ = (B[1]   - B[end-1])/2
        elseif i == 1
            dβ′ = (B[2]   - B[end]  )/2
        end

        ɛk² = sparse(1:N, 1:N, ɛ[:]*k², N, N)
        χk² = sparse(1:N, 1:N, D₀*γ(sim,k)*F.*F0[:]*k², N, N)
        β′∇ = sparse(1:N, 1:N, β′, N, N)*∇

        A = (∇²+ɛk²+χk²+2im*β′∇-β′²)\M
        A₀ += A*dβ′/complex(0,2π)
        A₁ += A*β′*dβ′/complex(0,2π)

    end

    P = svd(A₀)
    temp = findall(P.S .< rank_tol)
    if isempty(temp)
        @warn "need more evals. Rerun with larger nβ"
        return ComplexF64[NaN]
    else
        k = temp[1]-1
    end

    B = (P.U[:,1:k])'*A₁*(P.Vt[1:k,:])'*diagm(0 => 1 ./P.S[1:k])

    F = eigen(B)

    return F.values::Array{ComplexF64,1}
end
