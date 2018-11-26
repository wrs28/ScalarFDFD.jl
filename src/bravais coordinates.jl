################################################################################
# COORDINATE TRANSFORMATIONS
################################################################################
"""
    xb, yb = bravais_coordinates_unit_cell(x, y, domain_index, sys)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
sys.domains[domain_index].lattice
"""
function bravais_coordinates_unit_cell(x, y, domain::Int, system::System)
    return bravais_coordinates_unit_cell(x, y, system.domains[domain])
end
"""
    xb, yb = bravais_coordinates_unit_cell(x, y, domain)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
domain.lattice
"""
function bravais_coordinates_unit_cell(x, y, domain::Domain)
    return bravais_coordinates_unit_cell(x,y,domain.lattice)
end
"""
    xb, yb = bravais_coordinates_unit_cell(x, y, sim)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
sim.lat
"""
function bravais_coordinates_unit_cell(x, y, sim::Simulation)
    return bravais_coordinates_unit_cell(x,y,sim.lat)
end
"""
    xb, yb = bravais_coordinates_unit_cell(x, y, lattice)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
lattice
"""
function bravais_coordinates_unit_cell(x, y, lattice::Bravais)

    if (isinf(lattice.a) && iszero(lattice.β-π/2)) || (isinf(lattice.b) && iszero(lattice.α)) || (isinf(lattice.a) && isinf(lattice.b))
        return float(x), float(y)
    end

    p1, p2 = bravais_coordinates(x, y, lattice)
    p1 = mod(p1, lattice.a)
    p2 = mod(p2, lattice.b)

    xb = lattice.v1[1]*p1 + lattice.v2[1]*p2 + lattice.x0
    yb = lattice.v1[2]*p1 + lattice.v2[2]*p2 + lattice.y0

    return xb, yb
end


"""
    p1, p2 = bravais_coordinates(x, y, domain_index, sys)

coordinates in bravais frame specified by sys.domains[domain_index].lattice
(i.e. (x,y) = p1*v1 + p2*v2)
"""
function bravais_coordinates(x, y, domain::Int, system::System)
    return bravais_coordinates(x, y, system.domains[domain])
end
"""
    p1, p2 = bravais_coordinates(x, y, domain)

coordinates in bravais frame specified by domain.lattice
(i.e. (x,y) = p1*v1 + p2*v2)
"""
function bravais_coordinates(x, y, domain::Domain)
    bravais_coordinates(x, y, domain.lattice)
end
"""
    p1, p2 = bravais_coordinates(x, y, sim)

coordinates in bravais frame specified by sim.lat
(i.e. (x,y) = p1*v1 + p2*v2)
"""
function bravais_coordinates(x, y, sim::Simulation)
    return bravais_coordinates(x, y, sim.lat)
end
"""
    p1, p2 = bravais_coordinates(x,y,lattice)

coordinates in bravais frame (i.e. (x,y) = p1*v1 + p2*v2)
"""
function bravais_coordinates(x, y, lattice::Bravais)
    if (isinf(lattice.a) && iszero(lattice.β-π/2)) || (isinf(lattice.b) && iszero(lattice.α)) || (isinf(lattice.a) && isinf(lattice.b))
        return float(x-lattice.x0), float(y-lattice.y0)
    end

    x += -lattice.x0
    y += -lattice.y0

    rv1 = lattice.v1[1]*x + lattice.v1[2]*y
    rv2 = lattice.v2[1]*x + lattice.v2[2]*y

    p1 = (rv1 - rv2*lattice.cosθ)/lattice.sin²θ
    p2 = (rv2 - rv1*lattice.cosθ)/lattice.sin²θ
    return p1, p2
end
