#TODO: analyze_field for 1dim
#TODO: analyze_into_waveguides account for phase, for now it's not accounted for at all
#TODO: waveguide analysis only for perpendicular waveguides

"""
c = analyze_output(sim, k, ψ, m)
    s is the output coefficient in the mth channel

    S is constructed from s for unit input on each channel
"""
function analyze_field(sim::Simulation, k, ψ, m; direction=:out)
    if isempty(sim.sys.waveguides)
        c = analyze_into_angular_momentum(sim, k, ψ, m, direction)
    else
        c = analyze_into_waveguides(sim, k, ψ, m)
    end
    return c
end


"""
    analyze_into_angular_momentum(sim, k, ψ, m, direction)
"""
function analyze_into_angular_momentum(sim::Simulation, k, ψ, m, direction)

    nθ = NUMBER_OF_QUADRATURE_POINTS_ANG_MOM_ANALYSIS::Int
    θ = LinRange(0,2π,nθ)
    dθ = θ[2]-θ[1]

    # R is radius at which to interpolate
    R = (findmin(abs.(sim.bnd.∂Ω_tr)))[1] - findmax(sim.dis.dx)[1]
    X = R*cos.(θ)
    Y = R*sin.(θ)

    # interpolate wavefunction at r=R, result is P(θ)
    xs = sim.dis.x[1][1]:sim.dis.dx[1]:sim.dis.x[1][end]
    ys = sim.dis.x[2][1]:sim.dis.dx[2]:sim.dis.x[2][end]
    Ψ = CubicSplineInterpolation((xs,ys), reshape(ψ,sim.dis.N[1],sim.dis.N[2]))
    P = Ψ.(X,Y)

    q = sim.sct.channels[m].quantum_number

    if direction == :in
        c = quadrature(exp.(-1im*q*θ).*P, dθ)/(π*hankelh2(q,k*R))
    elseif direction == :out
        c = quadrature(exp.(-1im*q*θ).*P, dθ)/(π*hankelh1(q,k*R))
    else
        throw(ErrorException("Invalid direction $(direction). Must be one of :in or :out"))
    end
    return c
end


function analyze_into_waveguides(sim::Simulation, k, ψ, m)
    if any(ispc.(sim.sys.domains[get_waveguide_domains(sim,sim.sct.channels[m].waveguide)]))
        c = analyze_into_waveguides_pc(sim, k, ψ, m)
    else
        c = analyze_into_waveguides_planar(sim, k, ψ, m)
    end
    return c
end

"""
    c = analyze_into_waveguides_pc(sim, k, ψ, m)
"""
function analyze_into_waveguides_pc(sim::Simulation, k, ψ, m)

    xs = sim.dis.x[1][1] .+ (1:sim.dis.N[1])*sim.dis.dx[1]
    ys = sim.dis.x[2][1] .+ (1:sim.dis.N[2])*sim.dis.dx[2]
    utp = CubicSplineInterpolation((xs,ys),reshape(ψ[:,1],sim.dis.N[1],sim.dis.N[2]), bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    _, _, _, φ_itp, prop_const = incident_mode_pc(sim, k, m)

    wg_sim = extract_waveguide_simulation(sim,sim.sct.channels[m].waveguide)

    n_quad = 100

    if get_asymptote(sim, sim.sct.channels[m].waveguide) == :left
        xs = sim.bnd.∂Ω_tr[1,1] .+ LinRange(0, wg_sim.lat.a, n_quad)
        ys = sim.dis.x[2]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :right
        xs = sim.bnd.∂Ω_tr[1,2] .- LinRange(0, wg_sim.lat.a, n_quad)
        ys = sim.dis.x[2]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :bottom
        xs = sim.dis.x[1]
        ys = permutedims(sim.bnd.∂Ω_tr[2,1] .+ LinRange(0, wg_sim.lat.b, n_quad))
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :top
        xs = sim.dis.x[1]
        ys = permutedims(sim.bnd.∂Ω_tr[2,2] .- LinRange(0, wg_sim.lat.b, n_quad))
    else
        throw(ArgumentError("invalid asymptote $(get_asymptote(sim,sim.sct.channels[i].waveguide)) for channel $i"))
    end

    dx = xs[2]-xs[1]
    dy = ys[2]-ys[1]
    c = prop_const*quadrature(utp.(xs,ys).*conj(φ_itp.(xs,ys)), [dx,dy])::ComplexF64

    return c
end

"""
    analyze_into_waveguides(sim, k, ψ, m; direction)
"""
function analyze_into_waveguides_planar(sim::Simulation, k, ψ, m, direction; η_init=-3)

    side = sim.bnd.waveguides[sim.sct.channels[m].waveguide].side
    xs = sim.dis.x[1][1]:sim.dis.dx[1]:sim.dis.x[1][end]
    ys = sim.dis.x[2][1]:sim.dis.dx[2]:sim.dis.x[2][end]
    Ψ = CubicSplineInterpolation((xs,ys), reshape(ψ,sim.dis.N[1],sim.dis.N[2]))

    β, u_itp, collapsed_dim, full_dim = transverse_mode(sim, k, m)
    φ = u_itp(sim.dis.x[full_dim])
    if side == :l
        P = Ψ(sim.dis.x[1][1],sim.dis.x[2])
    elseif side == :r
        P = Ψ(sim.dis.x[1][end],sim.dis.x[2])
    elseif side == :b
        P = Ψ(sim.dis.x[1],sim.dis.x[2][1])
    elseif side == :t
        P = Ψ(sim.dis.x[1],sim.dis.x[2][end])
    else
        throw(ErrorException("Invalid side $(side)"))
    end

    c = sqrt(real(β))*quadrature((conj(φ).*P)[:], sim.dis.dx[full_dim])
    return c
end


"""
    flux, (left, right, bottom, top) = surface_flux(sim, ψ)
"""
function surface_flux(sim::Simulation,Ψ)

    flux = Array{Float64}(undef,size(Ψ,2))
    top = Array{Float64}(undef,size(flux))
    bottom = Array{Float64}(undef,size(flux))
    left = Array{Float64}(undef,size(flux))
    right = Array{Float64}(undef,size(flux))

    for i in 1:size(Ψ,2)
        ψ = reshape(Ψ[sim.dis.X_idx,i],sim.dis.N_tr[1],sim.dis.N_tr[2])
        ∇₁, ∇₂ = grad(sim.dis.N_tr,sim.dis.dx)
        ∇₁ψ = reshape(∇₁*ψ[:], sim.dis.N_tr[1]-1, sim.dis.N_tr[2])
        ∇₂ψ = reshape(∇₂*ψ[:], sim.dis.N_tr[1], sim.dis.N_tr[2]-1)
        kx = imag((∇₁ψ).*conj(ψ[1:end-1,:]+ψ[2:end,:])/2)
        ky = imag((∇₂ψ).*conj(ψ[:,1:end-1]+ψ[:,2:end])/2)
        bottom[i] = quadrature(ky[:,1],[sim.dis.dx[1]])
        top[i] = quadrature(ky[:,end],[sim.dis.dx[1]])
        left[i]  = quadrature(kx[1,:],[sim.dis.dx[2]])
        right[i]  = quadrature(kx[end,:],[sim.dis.dx[2]])
        flux[i] = top[i] + bottom[i] + right[i] + left[i]
    end

    return flux, (left,right,bottom,top)
end


"""
    absorption = compute_loss(sim, k, ψ)
"""
function bulk_absorption(sim::Simulation, k, ψ)

    absorption = Array{Float64}(undef, length(k))

    for i ∈ eachindex(k)

        Ψ = ψ[sim.dis.X_idx,i]
        A = imag((sim.sys.ε[sim.dis.X_idx]-1im*sim.tls.D₀*sim.sys.F[sim.dis.X_idx])*k[i]^2).*abs2.(Ψ)
        absorption[i] = quadrature(sim, A)
    end

    return absorption
end
