
"""
    lattice = Bravais(bnd::Boundary)

Rectangular lattice with same size as defined in `bnd` along periodic directions.
"""
function Bravais(bnd::Boundary)
    if :p ∈ bnd.bc[:,1] && bnd.bc[1,1] == bnd.bc[2,1]
        a = bnd.∂Ω[2,1]-bnd.∂Ω[1,1]
    elseif :p ∈ bnd.bc[:,1] && bnd.bc[1,1] !== bnd.bc[2,1]
        throw(ArgumentError("only one of bc[1,1] and bc[2,1] is :p"))
    else
        a = Inf
    end
    if :p ∈ bnd.bc[:,2] && bnd.bc[1,2] == bnd.bc[2,2]
        b = bnd.∂Ω[2,2]-bnd.∂Ω[1,2]
    elseif :p ∈ bnd.bc[:,2] && bnd.bc[1,2] !== bnd.bc[2,2]
        throw(ArgumentError("only one of bc[1,2] and bc[2,2] is :p"))
    else
        b = Inf
    end
    return Bravais(a=a, b=b)
end


"""
    lattice = Bravais(sim::Simulation)

Rectangular lattice with same size as defined in `sim.bnd` along periodic directions.
"""
function Bravais(sim::Simulation)
    return Bravais(sim.bnd)
end


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
