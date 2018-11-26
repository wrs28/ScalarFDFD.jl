################################################################################
## QUADRATURE
################################################################################
"""
    ∫f(z) dz = quadrature(sim, f; weight=:none)

multi-dimensional integral of f. really just sum over sites.
"""
function quadrature(sim::Simulation, f::AbstractArray{T,N}; weight=:none, k=1) where T where N

    if weight == :ε
        Σ = sum(f[:].*sim.sys.ε[:])*prod(sim.dis.dx)
    elseif weight == :ε_bl
        ε = ε_bl(sim; k=k)
        Σ = sum(f[:].*ε)*prod(sim.dis.dx)
    else
        Σ = sum(f[:])*prod(sim.dis.dx)
    end

    return Σ
end


function ε_bl(sim::Simulation; k=1)
    SA = absorbing_boundary_layers(sim, k)
    S1, S2 = pml_boundary_layers(sim, k)

    S1 = S2[1,1]*S2[2,1]
    S1 = sparse(1:sim.dis.N[1], 1:sim.dis.N[1], 1 ./S1.nzval)
    # S1 = (SA[1,1]*SA[2,1])^2*S1

    S2 = S2[1,2]*S2[2,2]
    S2 = sparse(1:sim.dis.N[2], 1:sim.dis.N[2], 1 ./S2.nzval)
    # S2 = (SA[1,2]*SA[2,2])^2*S2

    S = kron(S2,S1)
    SA = diag(kron(SA[1,2]+SA[2,2],sparse(1.0*I,sim.dis.N[1],sim.dis.N[1])) + kron(sparse(1.0*I,sim.dis.N[2],sim.dis.N[2]),SA[1,1]+SA[2,1]))

    return (sqrt.(S*sim.sys.ε[:]) + SA).^2
end
