################################################################################
### PLANAR WAVEGUIDES
################################################################################

# TODO: Test arbitrary profile waveguide, for now have only tested simple waveguide
# TODO: make waveguide addition per side, as with pc waveguides
"""
    waveguide_domains = planar_waveguide(regions::Regions; direction)
"""
function planar_waveguide(regions::Regions; direction)
    if direction ∈ [:x, :X, :l, :r, :lr, :L, :R, :LR, :RL, :rl, :horizontal, :h]
        bulk_domain_type = :bulk_waveguide_x
        which_asymptote1 = :left
        which_asymptote2 = :right
    elseif direction ∈ [:y, :Y, :u, :d, :ud, :D, :U, :UD, :DU, :du, :vertical, :v]
        bulk_domain_type = :bulk_waveguide_y
        which_asymptote1 = :bottom
        which_asymptote2 = :top
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end
    waveguide_number1 = 1#sort(unique([old_domains[i].which_waveguide for i ∈ eachindex(old_domains)]))[end]+1
    waveguide_number2 = 1#waveguide_number1+1
    wvg_bulk = build_domain(regions; is_in_domain=regions.functions[1], domain_params=regions.params[1], domain_type=bulk_domain_type)
    wvg1 = build_domain(regions; which_asymptote=which_asymptote1, which_waveguide=waveguide_number1)
    wvg2 = build_domain(regions; which_asymptote=which_asymptote2, which_waveguide=waveguide_number2)
    return vcat(wvg1, wvg2, wvg_bulk)
end
"""
    waveguide_domains = planar_waveguides(;width, index, direction, transverse_position)
"""
function planar_waveguide(;width, index, direction, transverse_position)
    if direction ∈ [:x, :X, :l, :r, :lr, :L, :R, :LR, :RL, :rl, :horizontal, :h]
        regions = add_rectangle(a=Inf, b=width, x0=-1e4, y0=transverse_position-width/2, θ=0, n₁=index)
    elseif direction ∈ [:y, :Y, :u, :d, :ud, :D, :U, :UD, :DU, :du, :vertical, :v]
        regions = add_rectangle(a=width, b=Inf, x0=transverse_position-width/2, y0=-1e4, θ=0, n₁=index)
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end
    return add_planar_waveguide(regions; direction=direction)
end


"""
    simulation = add_waveguide(sim; width, index, direction, position)
    simulation = add_waveguide(regions, sim; width, index, direction, position)

    system = add_waveguide(sys::System; width, index, direction, position)
    system = add_waveguide(regions::Regions, sys::System; width, index, direction, position)

    domains = add_waveguide(domains::Array{Domains}; width, index, direction, position)
    domains = add_waveguide(regions::Regions, domains::Array{Domains}; width, index, direction, position)
"""
function add_planar_waveguide(sim::Simulation; width, index, direction, position)
    sys = add_planar_waveguide(sim.sys; width=width, index=index, direction=direction, position=position)
    return Simulation(sim; :sys => sys)
end
function add_planar_waveguide(sys::System; width, index, direction, position)
    waveguide_domains = planar_waveguide(; width=width, index=index, direction=direction, position=position)
    return System(vcat(waveguide_domains,sys.domains))
end
function add_planar_waveguide(regions::Regions, sim::Simulation; direction)
    sys = add_planar_waveguide(regions, sim.sys; direction=direction)
    return Simulation(sim; :sys => sys)
end
function add_planar_waveguide(regions::Regions, sys::System; direction)
    waveguide_domains = planar_waveguide(regions; direction=direction)
    return System(vcat(waveguide_domains,sys.domains))
end


################################################################################
### PHOTONIC CRYSTAL WAVEGUIDES
################################################################################


#TODO: fix waveguide identifier
"""
    waveguide_domains = defect_waveguide(domain::Domain, x0, y0, direction, width)

domain should be of the :pc type.
"""
function defect_waveguide(domain::Domain; x0, y0, direction, width, number=1)
        xb0, yb0 = bravais_coordinates(x0, y0, domain)
        if direction ∈ [:u, :up, :n, :north, :t, :top]
            direction = :y
            N1_bulk = floor(Int,yb0/domain.lattice.b)
            N2_bulk = +1e4
            Nt = floor(Int,xb0/domain.lattice.a)
            N_width = ceil(Int,width/domain.lattice.a)
            bulk_domain_type = :bulk_pc_defect_waveguide_y
            which_asymptote = :top
        elseif direction ∈ [:b, :bottom, :s, :south, :d, :down]
            direction = :y
            N1_bulk = -1e4
            N2_bulk = floor(Int,yb0/domain.lattice.b)
            Nt = floor(Int,xb0/domain.lattice.a)
            N_width = ceil(Int,width/domain.lattice.a)
            bulk_domain_type = :bulk_pc_defect_waveguide_y
            which_asymptote = :bottom
        elseif direction ∈ [:r, :right, :e, :east]
            direction = :x
            N1_bulk = floor(Int,xb0/domain.lattice.a)
            N2_bulk = +1e4
            Nt = floor(Int,yb0/domain.lattice.b)
            N_width = ceil(Int,width/domain.lattice.b)
            bulk_domain_type = :bulk_pc_defect_waveguide_x
            which_asymptote = :right
        elseif direction ∈ [:l, :left, :w, :west]
            direction = :x
            N1_bulk = -1e4
            N2_bulk = floor(Int,xb0/domain.lattice.a)
            Nt = floor(Int,yb0/domain.lattice.b)
            N_width = ceil(Int,width/domain.lattice.b)
            bulk_domain_type = :bulk_pc_defect_waveguide_x
            which_asymptote = :left
        else
            throw(ArgumentError("invalid waveguide direction $(direction)"))
        end
        defect_bulk = add_line_defect(domain, direction, N1_bulk, N2_bulk, Nt; N_width=N_width)
        defect_bnd = add_line_defect(domain, direction, -1e4, 1e4, Nt; N_width=N_width)
        waveguide_number = number
        wvg_bulk = vcat(Domain(defect_bulk; :domain_type => bulk_domain_type))
        wvg = vcat(
                Domain(defect_bnd; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote),
                Domain(domain; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote)
                )

    return vcat(wvg, wvg_bulk)
end


"""
    sim = add_defect_waveguide(sim::Simulation; x0, y0, direction, width)
    sys = add_defect_waveguide(sys::System; x0, y0, direction, width)
"""
function add_defect_waveguide(sim::Simulation; x0, y0, direction, width)
    sys = add_defect_waveguide(sim.sys; x0=x0, y0=y0, direction=direction, width=width)
    return Simulation(sim; :sys => sys)
end
function add_defect_waveguide(sys::System; x0, y0, direction, width)
    domain = which_domain(x0,y0, Boundary([-Inf Inf;-Inf Inf],:d,:none,0), sys)
    waveguides = defect_waveguide(sys.domains[domain]; x0=x0, y0=y0, direction=direction, width=width)
    return System(vcat(waveguides,sys.domains))
end
