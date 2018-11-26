#TODO: consider use of Optim in CF rootfinding

################################################################################
# RESONANCE SOLVERS
################################################################################
"""
    k, ψ, η, conv = eig_knl(sim, k_init; F=[1], η_init=0, u_init=[], k_avoid=[0],
        disp_opt=false, tol=.5, max_count=15, max_iter=50, ka=0, kb=0)

eigenfrequency `k` with dispersion via root-finding of `η(k) - D₀γ(k)`.

`max_count`: the max number of CF states to compute at a time.
`max_iter`: maximum number of iterations to include in nonlinear solve
"""
function eig_knl(sim::Simulation, k_init::Number, η_init::Number, ka::Number, kb::Number, u_init::Array,
    k_avoid=[0], disp_opt::Bool=false, tol::Float64=.5, max_count::Int=15, max_iter::Int=50)

    ηt, ut = eig_cf(sim, k_init, 1, η_init, ka, kb, F, u_init)

    γp = sim.tls.γp; k₀ = sim.tls.k₀; D₀ = sim.tls.D₀; dx = sim.dis.dx
    function f!(fvec,z)

        k = complex(z[1],z[2])

        flag = true
        count = 1
        M = 1
        ind = Int

        while flag
            η_temp, u_temp = eig_cf(sim, k, M, ηt[1], ka, kb, ut[:,1])
            overlap = zeros(Float64,M)
            for i ∈ eachindex(overlap)
                overlap[i] = abs(quadrature(sim,ut[:,1].*sim.sys.F[:].*F.*u_temp[:,i]))
            end
            max_overlap, ind = findmax(overlap)
            if max_overlap > (1-tol)
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
            elseif count < max_count && max_overlap ≤ (1-tol)
                M += 1
            else
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
                @warn "overlap less than tolerance, skipped to neighboring TCF."
            end
            count += 1
        end
        fvec[1] = real((ηt[1]-D₀*γ(sim,k))/prod(k .- k_avoid))
        fvec[2] = imag((ηt[1]-D₀*γ(sim,k))/prod(k .- k_avoid))
    end

    z = nlsolve(f!, [real(float(k_init)),imag(float(k_init))]; iterations=max_iter, show_trace=disp_opt)
    k = complex(z.zero[1], z.zero[2])
    conv = converged(z)

    ηt[1] = D₀*γ(sim,k)
    η, ψ = eig_cf(sim, k, 1, ηt[1], ka, kb, F, ut[:,1])

    if !conv
        warn("no convergence for frequency $(k[1])")
    end
    return [k]::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}, η::Array{ComplexF64,1}, conv::Bool
end


"""
    k = eig_knl(sim, k, nk, radii; F=[1], r_min=.01, Nq=100, rank_tol=1e-8, ka=0, kb=0)

Eigenfrequency with dispersion, using contour integration. BC's and BL's  set
by `sim.bnd.bc`, `sim.bnd.bl`

Contour is centered on `k`, `radii` = [real radius, imag radius].

`nk` is an upper bound on the number of eigenfrequencies contained in the contour.

`Nq` is the number of contour quadrature points.
"""
function eig_knl(sim::Simulation, k::Number, nk::Int, radii::Union{Tuple,Array}, ka::Number, kb::Number, Nq::Int, rank_tol::Float64, r_min::Float64)

    ∇² = laplacian(sim, k; ka=ka, kb=kb)

    N = prod(sim.dis.N)
    A  = zeros(ComplexF64, N, nk)
    A₀ = zeros(ComplexF64, N, nk)
    A₁ = zeros(ComplexF64, N, nk)
    M = rand(N, nk)

    ϕ = 2π*range(0,length=Nq,step=1/Nq)
    Ω = k .+ complex.(radii[1]*cos.(ϕ), radii[2]*sin.(ϕ))

    θ =  angle(complex(sim.tls.k₀,-sim.tls.γp) - k)
    flag = abs(complex(sim.tls.k₀,-sim.tls.γp) - k) < rad.(radii[1],radii[2],θ)

    if flag
        AA  = zeros(ComplexF64, N, nk)
        AA₀ = zeros(ComplexF64, N, nk)
        AA₁ = zeros(ComplexF64, N, nk)
        r = 2r_min
        ΩΩ = complex(sim.tls.k₀,-sim.tls.γp) .+ complex.(rr*cos.(ϕ)/2,rr*sin.(ϕ)/2)
    end

    ε=sim.sys.ε; F0=sim.sys.F; D₀=sim.tls.D₀; dx = sim.dis.dx
    for i in 1:Nq

        k′ = Ω[i]
        k′² = k′^2

        if i > 1 && i < Nq
            dk′ = (Ω[i+1] - Ω[i-1]  )/2
        elseif i == Nq
            dk′ = (Ω[1]   - Ω[end-1])/2
        elseif i == 1
            dk′ = (Ω[2]   - Ω[end]  )/2
        end

        ɛk′² = sparse(1:N, 1:N, ɛ[:]*k′², N, N)
        χk′² = sparse(1:N, 1:N, D₀*γ(sim,k′).*F0[:]*k′², N, N)

        A = (∇²+ɛk′²+χk′²)\M
        A₀ += A*dk′/complex(0,2π)
        A₁ += A*k′*dk′/complex(0,2π)

        if flag
            kk′ = ΩΩ[i]
            kk′² = kk′^2
            if i > 1 && i < Nq
                dkk′ = (ΩΩ[i+1] - ΩΩ[i-1]  )/2
            elseif i == Nq
                dkk′ = (ΩΩ[1]   - ΩΩ[end-1])/2
            elseif i == 1
                dkk′ = (ΩΩ[2]   - ΩΩ[end]  )/2
            end
            ɛkk′² = sparse(1:N, 1:N, ɛ[:]*kk′², N, N)
            χkk′² = sparse(1:N, 1:N, D₀*γ(sim, kk′)*F.*F0[:]*kk′², N, N)

            AA = (∇²+ɛkk′²+χkk′²)\M
            AA₀ += AA*dkk′/complex(0,2π)
            AA₁ += AA*kk′*dkk′/complex(0,2π)
       end

    end

    if flag
        A₀ = A₀-AA₀
        A₁ = A₁-AA₁
    end

    P = svd(A₀)
    temp = findall(P.S .< rank_tol)
    if isempty(temp)
        @warn "need more evals. Rerun with larger nk"
        return ComplexF64[NaN]
    else
        k = temp[1]-1
    end

    B = (P.U[:,1:k])'*A₁*(P.Vt[1:k,:])'*diagm(0 => 1 ./P.S[1:k])

    F = eigen(B)

    return F.values::Array{ComplexF64,1}
end


################################################################################
#### AUXILLIARIES
################################################################################
"""
    rad(a,b,θ)

complex radius as a function of angle `θ` for an ellipse with semi-major axis `a`
and semi-minor axis `b`.

For use in contour eigensolver.
"""
function rad(a, b, θ)
    return b/sqrt(sin(θ)^2+(b/a)^2cos(θ)^2)
end
