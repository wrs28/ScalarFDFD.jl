"""
    k, ψ =  eig_kl(sim, k, nk=1, ka=0, kb=0, F=[1], ψ_init=[])
"""
function eig_kl(sim::Simulation, k::Number, nk::Int=1, ka=0, kb=0)

    ∇² = laplacian(sim, k; ka=ka, kb=kb)

    N = prod(sim.dis.N)
    ε = sim.sys.ε

    ɛ⁻¹ = sparse(1:N, 1:N, 1 ./ɛ[:], N, N)

    k², ψ, nconv, niter, nmult, resid = eigs(-ɛ⁻¹*∇², which = :LM, nev = nk, sigma = k^2)
    k = sqrt.(k²)

    for i = 1:nk
        𝒩² = quadrature(sim, abs2.(ψ[:,i]); weight=:ε_bl, k=k[i])
        ψ[:,i] = ψ[:,i]/sqrt(𝒩²)
    end

    return k::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}
end


"""
    η,u = eig_cf(sim, k, ncf=1, η_init, ka, kb, F, u_init)
"""
function eig_cf(sim::Simulation, k::Number, ncf::Int, η_init, ka, kb, F)

    k²= k^2

    ∇² = laplacian(sim, k; ka=ka, kb=kb)

    N = sim.dis.N; ε = sim.sys.ε; F0 = sim.sys.F;
    ɛk² = sparse(1:prod(N), 1:prod(N), ɛ[:]*k², prod(N), prod(N))
    sF  = sparse(1:prod(N), 1:prod(N), sign.(F.*F0[:] .+ MINIMUM_F::Float64), prod(N), prod(N))
    FF  = sparse(1:prod(N), 1:prod(N),  abs.(F.*F0[:] .+ MINIMUM_F::Float64), prod(N), prod(N))

    η, u, nconv, niter, nmult, resid = eigs(-sF*(∇²+ɛk²)./k², FF, which=:LM, nev=ncf, sigma=η_init)

    inds = sim.dis.X_idx
    if length(F) == 1
        F_temp = F.*F0[inds]
    else
        F_temp = F[inds].*F0[inds]
    end

    for ii = 1:ncf
        𝒩² = quadrature(sim, u[inds,ii].*F_temp.*u[inds,ii]; weight=:none)
        u[:,ii] = u[:,ii]/sqrt(𝒩²)
    end

    return η::Array{ComplexF64,1}, u::Array{ComplexF64,2}
end


################################################################################
### PLANAR WAVEGUIDE PROPAGATION CONSTANT SOLVER
################################################################################
"""
    k, ψ =  planar_kl(sim, β, k_init, nk=1; F=[1])
"""
function planar_kl(sim::Simulation, β::Number, k_init::Number, nk::Int=1)

    ∇² = ScalarFDFD.laplacian(sim, k_init)

    N = prod(sim.dis.N)
    ε = sim.sys.ε

    ɛ⁻¹ = sparse(1:N, 1:N, 1 ./ɛ[:], N, N)
    B² = sparse(complex(β^2,0)*I, N, N)

    k², ψ, nconv, niter, nmult, resid = eigs(-ɛ⁻¹*(∇²-B²), which = :LM, nev = nk, sigma = k_init^2)
    k = sqrt.(k²)

    for ii = 1:nk
        𝒩² = quadrature(sim, abs2.(ψ[inds,ii]))
        ψ[:,ii] = ψ[:,ii]/sqrt(𝒩²)
    end

    return k::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}
end
