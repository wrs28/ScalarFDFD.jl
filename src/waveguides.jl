# TODO: Test arbitrary profile waveguide, for now have only tested simple waveguide

################################################################################
### PLANAR WAVEGUIDES
################################################################################
"""
    new_sys = add_planar_waveguide(sys; width, index, direction, position)

    new_sim = add_planar_waveguide(sim; width, index, direction, position)
"""
function add_planar_waveguide(sys::System; width, index, direction, x0, y0, background_index=1)
    waveguide_number = length(sys.waveguides)+1
    waveguide_domains = planar_waveguide_domains(; width=width, index=index, direction=direction,
        x0=x0, y0=y0, waveguide_number=waveguide_number, background_index=background_index)
    return System(vcat(waveguide_domains,sys.domains))
end
function add_planar_waveguide(sim::Simulation; width, index, direction, x0, y0, background_index=sim.sys.domains[end].n₁)
    sys = add_planar_waveguide(sim.sys; width=width, index=index, direction=direction, x0=x0, y0=y0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end

"""
    new_sim = add_planar_waveguides(sim, dim; width, index, position)
"""
function add_planar_waveguides(sim::Simulation, dim::Int; width, index, position, background_index=sim.sys.domains[end].n₁)
    if dim == 1
        return add_horizontal_planar_waveguides(sim; width=width, index=index, y0=position, background_index=background_index)
    elseif dim == 2
        return add_vertical_planar_waveguides(sim; width=width, index=index, x0=position, background_index=background_index)
    else
        throw(ArgumentError("invalid dimension $dim"))
    end
end

"""
    new_sim = add_horizontal_planar_waveguides(sim; width, index, y0)

    new_sys = add_horizontal_planar_waveguides(sys; width, index, y0)
"""
function add_horizontal_planar_waveguides(sys::System; width, index, y0, background_index=1)
    sys = add_planar_waveguide(sys; width=width, index=index, direction=:h, x0=0, y0=y0, background_index=background_index)
    return sys
end
function add_horizontal_planar_waveguides(sim::Simulation; width, index, y0, background_index=sim.sys.domains[end].n₁)
    sys = add_horizontal_planar_waveguides(sim.sys, width=width, index=index, y0=y0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end


"""
    new_sim = add_vertical_planar_waveguides(sim; width, index, x0)

    new_sys = add_vertical_planar_waveguides(sys; width, index, x0)
"""
function add_vertical_planar_waveguides(sys::System; width, index, x0, background_index=1)
    sys = add_planar_waveguide(sys; width=width, index=index, direction=:v, x0=x0, y0=0, background_index=background_index)
    return sys
end
function add_vertical_planar_waveguides(sim::Simulation; width, index, x0, background_index=sim.sys.domains[end].n₁)
    sys = add_vertical_planar_waveguides(sim.sys; width=width, index=index, x0=x0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end


"""
    waveguide_domains = planar_waveguides(;width, index, direction, x0, y0, waveguide_number)
"""
function planar_waveguide_domains(;width, index, direction, x0=0, y0=0, waveguide_number, background_index=1)
    if direction ∈ [:ew, :lr, :h, :horizontal, :x]
        bulk_domain_type = :bulk_planar_waveguide_x
        which_asymptote1 = :left
        which_asymptote2 = :right
        a = 1e5
        domain = rectangle_domain(a=a, b=width, x0=x0-a/2, y0=y0-width/2, n₁=index)
    elseif direction ∈ [:ns, :ud, :tb, :v, :vertical, :y]
        bulk_domain_type = :bulk_planar_waveguide_y
        which_asymptote1 = :top
        which_asymptote2 = :bottom
        b=1e5
        domain = rectangle_domain(a=width, b=b, x0=x0-width/2, y0=y0-b/2, n₁=index)
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end

    bgd1 = Domain(; :domain_type => :planar_waveguide_background, :n₁ => background_index, :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1)
    bgd2 = Domain(; :domain_type => :planar_waveguide_background, :n₁ => background_index, :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2)
    wvg1 = Domain(domain; :domain_type => :planar_waveguide, :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1)
    wvg2 = Domain(domain; :domain_type => :planar_waveguide, :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2)
    wvg_bulk = Domain(domain; :domain_type => bulk_domain_type)

    return vcat(bgd1, bgd2, wvg1, wvg2, wvg_bulk)
end


"""
    new_sim = add_planar_waveguide(sim, side::Symbol; index=1)

    new_sys = add_planar_waveguide(sys, side::Symbol; index=1)
"""
function add_planar_waveguide(sim, side::Symbol; index=1)
    throw(ErrorException("this is intended for metallic waveguides, haven't done yet"))
    if side == :top
        direction = :n
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[2,2]
        width = sim.bnd.∂Ω[2,1]-sim.bnd.∂Ω[1,1]
    elseif side == :bottom
        direction = :s
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,1]-sim.bnd.∂Ω[1,1]
    elseif side == :left
        direction = :w
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,2]-sim.bnd.∂Ω[1,2]
    elseif side == :right
        direction = :e
        x0 = sim.bnd.∂Ω[2,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,2]-sim.bnd.∂Ω[1,2]
    else
        throw(ArgumentError("invalide side $side, must be one of :top, :bottom, :left, :right"))
    end
    return add_planar_waveguide(simsys; width=width, index=index, direction=direction, x0=x0, y0=y0)
end

################################################################################
### HALF-SPACE WAVEGUIDES (for 1-dim)
################################################################################
"""
    waveguide_domains = planar_waveguides(;index, direction, position, waveguide_number)
"""
function halfspace_waveguide_domains(;index, direction, position, waveguide_number)
    if direction ∈ [:w, :l, :left, :west]
        which_asymptote = :left
        lattice = Bravais(a=0)
    elseif direction ∈ [:e, :r, :right, :east]
        which_asymptote = :right
        lattice = Bravais(a=0)
    elseif direction ∈ [:n, :u, :t, :north, :up, :top]
        which_asymptote = :top
        lattice = Bravais(b=0)
    elseif direction ∈ [:s, :b, :d, :south, :bottom, :down]
        which_asymptote = :bottom
        lattice = Bravais(b=0)
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end

    wvg = build_domain(; which_asymptote=which_asymptote, which_waveguide=waveguide_number, domain_type=:halfspace_waveguide, n₁ = index, lattice=lattice)

    return wvg
end


"""
    sys = add_halfspace_waveguide(sys; index, direction, position)
"""
function add_halfspace_waveguide(sys::System; index, direction, position, waveguide_number=length(sys.waveguides)+1)
    waveguide_domains = halfspace_waveguide_domains(; index=index, direction=direction, position=position, waveguide_number=waveguide_number)
    return System(vcat(waveguide_domains,sys.domains))
end
"""
    sim = add_halfspace_waveguide(sim; index, direction)
"""
function add_halfspace_waveguide(sim::Simulation; index, direction)
    waveguide_number=length(sim.sys.waveguides)+1
    if sim.dis.N[1]==1
        position = direction ∈ [:n, :north, :t, :top] ? sim.bnd.∂Ω_tr[2,2] : sim.bnd.∂Ω_tr[1,2]
    elseif sim.dis.N[2]==1
        position = direction ∈ [:e, :east, :r, :right] ? sim.bnd.∂Ω_tr[2,1] : sim.bnd.∂Ω_tr[1,1]
    else
        throw(ArgumentError("halfspace waveguides only for one-dimensional-systems"))
    end
    sys = add_halfspace_waveguide(sim.sys; index=index, direction=direction, position=position, waveguide_number=waveguide_number)
    return Simulation(sim; :sys => sys)
end

"""
    sim = add_halfspace_waveguide(sim, side::Symbol; width, index)
"""
function add_halfspace_waveguide(sim::Simulation, side::Symbol; index)
    if side == :bottom
        direction = :s
    elseif side == :top
        direction = :n
    elseif side == :left
        direction = :w
    elseif side == :right
        direction = :e
    else
        throw(ArgumentError("invalid side $side, must be one of :bottom, :top, :left, :right"))
    end
    return add_halfspace_waveguide(sim; index=index, direction=direction)
end
"""
    sim = add_halfspace_waveguides(sim; index=1)
"""
function add_halfspace_waveguides(sim::Simulation; index=1)
    if sim.dis.N[1]==1
        sim = add_halfspace_waveguide(sim, :bottom; index=index)
        sim = add_halfspace_waveguide(sim, :top; index=index)
    elseif sim.dis.N[2]==1
        sim = add_halfspace_waveguide_left(sim, :left; index=index)
        sim = add_halfspace_waveguide_right(sim, :right; index=index)
    else
        throw(ArgumentError("halfspace waveguide only for 1-dim"))
    end
    return sim
end


################################################################################
### PHOTONIC CRYSTAL WAVEGUIDES
################################################################################
"""
    waveguide_domains = defect_waveguide(domain::Domain, x0, y0, direction, width)

domain should be of the :pc type.
"""
function pc_waveguide_domains(pc_domain; width, direction, x0=0, y0=0, waveguide_number)

        xb0, yb0 = ScalarFDFD.bravais_coordinates(x0, y0, pc_domain)
        if direction ∈ [:u, :up, :n, :north, :t, :top]
            axis = :y
            N1_bulk = floor(Int,yb0/pc_domain.lattice.b)
            N2_bulk = +1e4
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote = :top
        elseif direction ∈ [:b, :bottom, :s, :south, :d, :down]
            axis = :y
            N1_bulk = -1e4
            N2_bulk = floor(Int,yb0/pc_domain.lattice.b)
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote = :bottom
        elseif direction ∈ [:r, :right, :e, :east]
            axis = :x
            N1_bulk = floor(Int,xb0/pc_domain.lattice.a)
            N2_bulk = +1e4
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote = :right
        elseif direction ∈ [:l, :left, :w, :west]
            axis = :x
            N1_bulk = -1e4
            N2_bulk = floor(Int,xb0/pc_domain.lattice.a)
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote = :left
        elseif direction ∈ [:h, :horizontal, :x]
            axis = :x
            N1_bulk = -1e4
            N2_bulk = +1e4
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote1 = :left
            which_asymptote2 = :right
        elseif direction ∈ [:v, :vertical, :y]
            axis = :y
            N1_bulk = -1e4
            N2_bulk = +1e4
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote1 = :bottom
            which_asymptote2 = :top
        else
            throw(ArgumentError("invalid waveguide direction $(direction)"))
        end

        defect_bulk = line_defect_domain(pc_domain, axis, N1_bulk, N2_bulk, Nt; N_width=N_width)
        wvg_bulk = vcat(Domain(defect_bulk; :domain_type => bulk_domain_type))

        defect_bnd = line_defect_domain(pc_domain, axis, -1e4, 1e4, Nt; N_width=N_width)

        if direction ∈ [:v, :vertical, :y, :x, :horizontal, :h]
            wvg = vcat(
                    Domain(defect_bnd; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1, :domain_type => :pc_waveguide_background),
                    Domain(defect_bnd; :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2, :domain_type => :pc_waveguide_background)
                    )
        else
            wvg = vcat(
                    Domain(defect_bnd; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote, :domain_type => :pc_waveguide_background)
                    )

        end

    return vcat(wvg, wvg_bulk)
end


"""
    sim = add_defect_waveguide(sim::Simulation; x0, y0, direction, width)
    sys = add_defect_waveguide(sys::System; x0, y0, direction, width)
"""
function add_pc_waveguide(sim::Simulation; x0=0, y0=0, direction, width)
    sys = add_defect_waveguide(sim.sys; x0=x0, y0=y0, direction=direction, width=width)
    return Simulation(sim; :sys => sys)
end
function add_pc_waveguide(sys::System; x0=0, y0=0, direction, width)
    waveguide_number = length(sys.waveguides)+1
    temp_sys = System(sys.domains[.!ScalarFDFD.iswaveguide.([sys.domains[i] for i ∈ eachindex(sys.domains)]) .& .!ScalarFDFD.isbulkwaveguide.([sys.domains[i] for i ∈ eachindex(sys.domains)])])
    domain = ScalarFDFD.which_domain(x0,y0, Boundary([-Inf Inf;-Inf Inf],:d,:none,0), temp_sys)
    waveguide_domains = ScalarFDFD.pc_waveguide_domains(temp_sys.domains[domain]; x0=x0, y0=y0, direction=direction, width=width, waveguide_number=waveguide_number)
    return System(vcat(waveguide_domains,sys.domains))
end
