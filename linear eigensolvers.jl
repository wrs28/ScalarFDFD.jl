#TODO: revisit propagation constant solver
"""
    k, ψ =  eig_kl(sim, k, nk=1; F=[1], ψ_init=[], ka=0, kb=0)
"""
function eig_kl(sim::Simulation, k::Number, nk::Int=1, ka=0, kb=0, F=[1], ψ_init=[])

    ∇² = laplacian(sim, k; ka=ka, kb=kb)

    N = prod(sim.dis.N)
    ε = sim.sys.ε
    F0 = sim.sys.F
    D₀ = sim.tls.D₀

    ɛ⁻¹ = sparse(1:N, 1:N, 1 ./(ɛ[:] .- 1im*D₀.*F.*F0[:]), N, N)

    if isempty(ψ_init)
        k², ψ, nconv, niter, nmult, resid = eigs(-ɛ⁻¹*∇², which = :LM, nev = nk, sigma = k^2)
    else
        k², ψ, nconv, niter, nmult, resid = eigs(-ɛ⁻¹*∇², which = :LM, nev = nk, sigma = k^2, v0 = ψ_init)
    end

    inds = sim.dis.X_idx
    if length(F) == 1
        ε_temp = ɛ[inds]-1im*D₀.*F.*F0[inds]
    else
        ε_temp = ɛ[inds]-1im*D₀.*F[inds].*F0[inds]
    end

    for ii = 1:nk
        ψ[:,ii] = ψ[:,ii]/sqrt(quadrature(sim, ε_temp.*abs2.(ψ[inds,ii])))
    end

    return sqrt.(k²)::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}
end


"""
    η,u = eig_cf(sim, k, ncf=1; F=[1], η_init=0, u_init=[], ka=0, kb=0)
"""
function eig_cf(sim::Simulation, k::Number, ncf::Int=1, η_init=0, ka=0, kb=0, F=[1], u_init=[])

    k²= k^2

    ∇² = laplacian(sim, k; ka=ka, kb=kb)

    N = sim.dis.N; ε = sim.sys.ε; F0 = sim.sys.F
    ɛk² = sparse(1:prod(N), 1:prod(N), ɛ[:]*k², prod(N), prod(N))
    sF  = sparse(1:prod(N), 1:prod(N), sign.(F.*F0[:] .+ MINIMUM_F::Float64), prod(N), prod(N))
    FF  = sparse(1:prod(N), 1:prod(N),  abs.(F.*F0[:] .+ MINIMUM_F::Float64), prod(N), prod(N))

    if isempty(u_init)
        η, u, nconv, niter, nmult, resid = eigs(-sF*(∇²+ɛk²)./k², FF, which=:LM, nev=ncf, sigma=η_init)
    elseif !isempty(u_init)
        η, u, nconv, niter, nmult, resid = eigs(-sF*(∇²+ɛk²)./k², FF, which=:LM, nev=ncf, sigma=η_init, v0=u_init)
    end

    inds = sim.dis.X_idx
    if length(F) == 1
        F_temp = F.*F0[inds]
    else
        F_temp = F[inds].*F0[inds]
    end

    for ii = 1:ncf
        u[:,ii] = u[:,ii]/sqrt(quadrature(sim, u[inds,ii].*F_temp.*u[inds,ii]))
    end

    return η::Array{ComplexF64,1}, u::Array{ComplexF64,2}
end
