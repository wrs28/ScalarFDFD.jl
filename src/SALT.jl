#TODO: other SALT specific routines here

"""
    γ(sim, k) is the lorentzian gain curve
"""
function γ(sim::Simulation, k)::ComplexF64
    return sim.tls.γp/(k-sim.tls.k₀+1im*sim.tls.γp)
end
