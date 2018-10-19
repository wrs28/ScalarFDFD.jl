"""
    1. domain = build_pc(lattice, domain)
    2. domain = build_pc(lattice, regions=Regions(); is_in_domain=whole_domain, domain_params=[], n₁=1, n₂=0, F=0)

1. build a photonic crystal out of a domain by tiling whatever part of the domain
falls within the cell defined by `lattice`

2. similar to (1) but can build out of a Regions object instead.
"""
function build_pc(lattice::Bravais, domain::Domain)
    return build_pc(lattice, Regions=Regions(domain);
                    is_in_domain=domain.is_in_domain, domain_params=domain.domain_params,
                    n₁=domain.n₁, n₂=domain.n₂, F=domain.F)
end
function build_pc(lattice::Bravais, regions::Regions=Regions();
                is_in_domain=whole_domain, domain_params=[], n₁=1, n₂=0, F=0)
    return build_domain(regions; is_in_domain=is_in_domain, domain_params=domain_params,
                            n₁=n₁, n₂=n₂, F=F, domain_type=:pc, lattice=lattice)
end



"""
    domain = add_circle(domain, reference_point: R, x0=0, y0=0, n₁=1, n₂=0, F=0)

add a circular region to photonic crystal domain, referenced to `reference_point`,
which can be one of :center, :corner or :corners, :bottom_edge or :bottom or :top_edge
or :top, :left_edge or :left or :right_edge or :right, :southwest, :south, :southeast,
:east, :northeast, :north, :northwest, :west

this creates a circle which wraps around the edges of the cell, so that the final
structure has circles, not pie-wedges
"""
function add_circle(domain::Domain, reference::Symbol; R, x0=0, y0=0, n₁=1, n₂=0, F=0)
    v1 = domain.lattice.v1
    v2 = domain.lattice.v2
    a = domain.lattice.a
    b = domain.lattice.b

    if reference == :center
        x0_new = x0 - v1[1]*a/2 + v2[1]*b/2
        y0_new = y0 - v1[2]*a/2 + v2[2]*b/2
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = x0 + v1[1]*a/2 - v2[1]*b/2
        y0_new = y0 + v1[2]*a/2 - v2[2]*b/2
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = x0 + 3v1[1]*a/2 + v2[1]*b/2
        y0_new = y0 + 3v1[2]*a/2 + v2[2]*b/2
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = x0 + v1[1]*a/2 + 3v2[1]*b/2
        y0_new = y0 + v1[2]*a/2 + 3v2[2]*b/2
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = x0 + v1[1]*a/2 + v2[1]*b/2
        y0_new = y0 + v1[2]*a/2 + v2[2]*b/2
    elseif reference ∈ [:corner, :corners]
        x0_new = v1[1]*a
        y0_new = v1[2]*a
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = v2[1]*b
        y0_new = v2[2]*b
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = v1[1]*a + v2[1]*b
        y0_new = v1[2]*a + v2[2]*b
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = 0
        y0_new=0
    elseif reference ∈ [:bottom_edge, :top_edge, :bottom, :top]
        x0_new = v1[1]*a/2 + v2[1]*b
        y0_new = v1[2]*a/2 + v2[2]*b
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = v1[1]*a/2
        y0_new = v1[2]*a/2
    elseif reference ∈ [:left_edge, :right_edge, :left, :right]
        x0_new = v1[1]*a + v2[1]*b/2
        y0_new = v1[2]*a + v2[2]*b/2
        domain = add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)

        x0_new = v2[1]*b/2
        y0_new = v2[2]*b/2
    elseif reference == :southwest
        x0_new = x0
        y0_new = y0
    elseif reference == :south
        x0_new = x0 + v1[1]*a/2
        y0_new = y0 + v1[1]*a/2
    elseif reference == :southeast
        x0_new = x0 + v1[1]*a
        y0_new = y0 + v1[2]*a
    elseif reference ∈ :east
        x0_new = x0 + v1[1]*a + v2[1]*b/2
        y0_new = y0 + v1[2]*a + v2[2]*b/2
    elseif reference == :northeast
        x0_new = x0 + v1[1]*a + v2[1]*b
        y0_new = y0 + v1[2]*a + v2[2]*b
    elseif reference == :north
        x0_new = x0 + v1[1]*a/2 + v2[1]*b
        y0_new = y0 + v1[2]*a/2 + v2[2]*b
    elseif reference == :northwest
        x0_new = x0 + v2[1]*b
        y0_new = y0 + v2[2]*b
    elseif reference == :west
        x0_new = x0 + v2[1]*b/2
        y0_new = y0 + v2[2]*b/2
    else
        throw(ArgumentError("Invalid unit cell position reference $(reference).
        Should be one of :center, :corner or a cardinal direction"))
    end
    return add_circle(domain; R=R, x0=x0_new, y0=y0_new, n₁=n₁, n₂=n₂, F=F)
end

"""
    1. defect_domain = add_defect(domain, Na, Nb, regions=Region())
    2. new_system = add_defect(sys, domain_number, Na, Nb, regions=Region())

1. `defect_domain` is a domain that will replace the (Na,Nb) cell of the photonic
crystal `domain` (referenced to the `x0`, `y0` fields of `dom.lattice`)

2. similar to (1), but adds a defect to `sys.domains[domain_number]`

in both cases the interior of the defect is built from `regions`

"""
function add_defect(domain::Domain, Na::Int, Nb::Int, regions::Regions=Regions())
    defect_domain = add_line_defect(domain, :x, Na, Na, Nb, regions)
    return defect_domain
end
function add_defect(domain::Domain, Na::Int, Nb::Array{Int}, regions::Regions=Regions())
    return add_defect.(Ref(domain), Na, Nb, Ref(regions))
end
function add_defect(domain::Domain, Na::Array{Int}, Nb, regions::Regions=Regions())
    return add_defect.(Ref(domain), Na, Nb, Ref(regions))
end
function add_defect(sys::System, domain_num, Na, Nb, regions::Regions=Regions())
    domains = sys.domains
    defect_domain = add_defect(domains[domain_num], Na, Nb, regions)
    return System([defect_domain, domains])
end


# TODO: make defect sites non-periodic
"""
    defect_domain = add_line_defect(domain, direction, N_start, N_stop, N_transverse, regions=Regions()); N_step=1, N_width=1)

adds a line of defect sites, starting at site number `(N_start, N_transverse)` or
`(N_transverse, N_start)`, depending on the `direction`.

What happens inside can be specified by `regions` (for the time being `regions` is
repeated inside each defect cell).

Optional arguments:

`N_step` enables skipping of sites to add a defect to

`N_width` sets the width (in the transverse direction) of the line defect
"""
function add_line_defect(domain::Domain, direction, N_start, N_stop, N_transverse, regions::Regions=Regions(); N_step=1, N_width=1)
    if direction == :x
        Na_start = N_start
        Na_step = N_step
        Na_stop = N_stop
        Nb_start = N_transverse
        Nb_step = 1
        Nb_stop = N_transverse+N_width-1
    elseif direction == :y
        Nb_start = N_start
        Nb_step = N_step
        Nb_stop = N_stop
        Na_start = N_transverse
        Na_step = 1
        Na_stop = N_transverse+N_width-1
    else
        throw(ArgumentError("Invalid direction $(direction)."))
    end
    lattice = domain.lattice
    n₁ = domain.n₁
    n₂ = domain.n₂
    F = domain.F

    if Na_start==Na_stop && Nb_start==Nb_stop
        defect_type = :point_defect
    else
        if N_step > 1
            defect_type = :periodic_line_defect
        else
            defect_type = :simple_line_defect
        end
    end

    regions = deepcopy(regions)
    for i ∈ eachindex(regions.functions)
        regions.params[i][1] += Na_start*lattice.a*lattice.v1[1] * Nb_start*lattice.b*lattice.v2[1]
        regions.params[i][2] += Na_start*lattice.a*lattice.v1[2] * Nb_start*lattice.b*lattice.v2[2]
    end
    defect_domain = build_domain(regions; is_in_domain=is_in_cell,
                        domain_params=[Na_start, Na_stop, Na_step, Nb_start, Nb_stop, Nb_step],
                        n₁=n₁, n₂=n₂, F=F, domain_type=defect_type, lattice=lattice)
    return defect_domain
end


"""
    bool = is_in_cell(x, y, domain_number, bnd, sys)

is the point (`x`,`y`) in `sys.domains[domain_number]`?

`bnd` is just along for the ride
"""
function is_in_cell(x, y, idx, bnd::Boundary, sys::System)
    Na_start = Int(sys.domains[idx].domain_params[1])
    Na_stop = Int(sys.domains[idx].domain_params[2])
    Na_step = Int(sys.domains[idx].domain_params[3])
    Nb_start = Int(sys.domains[idx].domain_params[4])
    Nb_stop = Int(sys.domains[idx].domain_params[5])
    Nb_step = Int(sys.domains[idx].domain_params[6])

    if isinf(sys.domains[idx].lattice.a) && isinf(sys.domains[idx].lattice.b)
        return true
    end

    x += -sys.domains[idx].lattice.x0
    y += -sys.domains[idx].lattice.y0

    p1, p2 = bravais_coordinates(x, y, sys.domains[idx].lattice)

    pa = floor(Int,p1/sys.domains[idx].lattice.a)
    pb = floor(Int,p2/sys.domains[idx].lattice.b)

    return pa ∈ Na_start:Na_step:Na_stop && pb ∈ Nb_start:Nb_step:Nb_stop
end
