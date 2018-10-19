# TODO: parallel solves for nonlinear smatrix and eigensolvers (except for contour, it's done)
#TODO: fix parametric solver for cf, check for k
#TODO: documentation for CF solver (fix)

################################################################################
### EIGENVALUE SOLVERS OVER PARAMETERS, WRAPPERS
################################################################################
"""
    K = eig_k(sim, k, nk, systems::Array{System}, twolevelsystems::Array{TwoLevelSystem},
    lattices::Array{Bravais}; ka=0, kb=0, F=[1], disp_opt=true)

computes spectra of parametrically tuned simulations in parallel.

`systems` is an array of `System` objects, each of which will become an input to
a simulation. Same for `twolevelsystems` and `lattices`.

These three arrays can appear in any order (so long as they appear after `nk`),
and any two can be unspecified (e.g. if you are just tuning index but not lattice
or pump).

`K` has the same shape as `systems`, and the three arrays must all have the same
shape.
"""
function eig_k(sim::Simulation, k::Number, nk::Int,
    systems::Array{System,M},
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(systems)],
    lattices::Array{Bravais,M}=[sim.lat for i ∈ CartesianIndices(systems)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_p(eig_klp, sim, k, nk, systems, twolevelsystems, lattices, ka, kb, 0, F, ψ_init, disp_opt)
    return K
end
function eig_k(sim::Simulation, k::Number, nk::Int,
    systems::Array{System,M},
    lattices::Array{Bravais,M},
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(systems)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_k(sim, k, nk, systems, twolevelsystems, lattices; ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_k(sim::Simulation, k::Number, nk::Int,
    twolevelsystems::Array{TwoLevelSystem,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(twolevelsystems)],
    lattices::Array{Bravais,M}=[sim.lat for i ∈ CartesianIndices(twolevelsystems)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_k(sim, k, nk, systems, twolevelsystems, lattices; ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_k(sim::Simulation, k::Number, nk::Int,
    twolevelsystems::Array{TwoLevelSystem,M},
    lattices::Array{Bravais,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(twolevelsystems)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_k(sim, k, nk, systems, twolevelsystems, lattices; ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_k(sim::Simulation, k::Number, nk::Int,
    lattices::Array{Bravais,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(lattices)],
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(lattices)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_k(sim, k, nk, systems, twolevelsystems, lattices; ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_k(sim::Simulation, k::Number, nk::Int,
    lattices::Array{Bravais,M},
    twolevelsystems::Array{TwoLevelSystem,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(lattices)];
    ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_k(sim, k, nk, systems, twolevelsystems, lattices; ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end


"""
    K = eig_cf(sim, k, nk, systems::Array{System}, twolevelsystems::Array{TwoLevelSystem},
    lattices::Array{Bravais}; ka=0, kb=0, F=[1], disp_opt=true)

computes spectra of parametrically tuned simulations in parallel.

`systems` is an array of `System` objects, each of which will become an input to
a simulation. Same for `twolevelsystems` and `lattices`.

These three arrays can appear in any order (so long as they appear after `nk`),
and any two can be unspecified (e.g. if you are just tuning index but not lattice
or pump).

`K` has the same shape as `systems`, and the three arrays must all have the same
shape.
"""
function eig_cf(sim::Simulation, k, ncf::Int,
    systems::Array{System,M},
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(systems)],
    lattices::Array{Bravais,M}=[sim.lat for i ∈ CartesianIndices(systems)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_p(eig_cfp, sim, k, ncf, systems, twolevelsystems, lattices, η_init, ka, kb, 0, F, ψ_init, disp_opt)
    return K
end
function eig_cf(sim::Simulation, k, ncf::Int,
    systems::Array{System,M},
    lattices::Array{Bravais,M},
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(systems)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_cf(sim, k, ncf, systems, twolevelsystems, lattices; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_cf(sim::Simulation, k, ncf::Int,
    twolevelsystems::Array{TwoLevelSystem,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(twolevelsystems)],
    lattices::Array{Bravais,M}=[sim.lat for i ∈ CartesianIndices(twolevelsystems)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_cf(sim, k, ncf, systems, twolevelsystems, lattices; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_cf(sim::Simulation, k, ncf::Int,
    twolevelsystems::Array{TwoLevelSystem,M},
    lattices::Array{Bravais,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(twolevelsystems)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_cf(sim, k, ncf, systems, twolevelsystems, lattices; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_cf(sim::Simulation, k, ncf::Int,
    lattices::Array{Bravais,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(lattices)],
    twolevelsystems::Array{TwoLevelSystem,M}=[sim.tls for i ∈ CartesianIndices(lattices)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_cf(sim, k, ncf, systems, twolevelsystems, lattices; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end
function eig_cf(sim::Simulation, k, ncf::Int,
    lattices::Array{Bravais,M},
    twolevelsystems::Array{TwoLevelSystem,M},
    systems::Array{System,M}=[sim.sys for i ∈ CartesianIndices(lattices)];
    η_init=0, ka=0, kb=0, F=[1], ψ_init=[], disp_opt=true) where M

    K = eig_cf(sim, k, ncf, systems, twolevelsystems, lattices; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init, disp_opt=disp_opt)
    return K
end


################################################################################
### EIGENVALUE SOLVERS OVER PARAMETERS, CORE
################################################################################
"""
    E = eig_p(eig_fun, sim, k, systems, twolevelsystems, lattices, ka, kb, η_init, F, ψ_init, disp_opt)
"""
function eig_p(eig_fun::Function, sim::Simulation, k, nk::Int,
    systems::Array{System,M},
    twolevelsystems::Array{TwoLevelSystem,M},
    lattices::Array{Bravais,M},
    ka=0, kb=0, η_init=0, F=[1], ψ_init=[], disp_opt=true) where M

    E = Array{ComplexF64}(undef, nk, size(systems)...)

    sim = Simulation(sim; sys=systems[1], tls=twolevelsystems[1], lattice=lattices[1])
    if eig_fun == eig_klp
        E[:,fill(1,M)...] = eig_k(sim, k, nk; ka=ka, kb=kb, F=F, ψ_init=ψ_init)[1]
    elseif eig_fun == eig_cfp
        E[:,fill(1,M)...] = eig_cf(sim, k[1], nk; η_init=η_init, ka=ka, kb=kb, F=F, ψ_init=ψ_init)[1]
    else
        throw(ArgumentError("parametric parallization only available for linear eig_k and eig_cf at present."))
    end

    for e_dim ∈ 1:M
        e_inds = Array{AbstractArray}(undef,M+1)
        e_inds[1] = 1:nk
        pres_inds = Array{AbstractArray}(undef,M)
        for i ∈ 1:M
            if i < e_dim
                pres_inds[i] = 1:size(systems, i)
                e_inds[i+1] = 1:size(systems,i)
            elseif i == e_dim
                pres_inds[i] = 1:size(systems, i)
                e_inds[i+1] = 1:1
            else
                pres_inds[i] = 1:1
                e_inds[i+1] = 1:1
            end
        end
        pres_inds = (pres_inds...,)
        e_inds = (e_inds...,)

        sims = [sim for i ∈ CartesianIndices(E[e_inds...])]
        Es = [E[e_inds...][i] for i ∈ CartesianIndices(E[e_inds...])]
        syss = [systems[pres_inds...][Tuple(i)[2:e_dim]...,pres_inds[e_dim:end]...] for i ∈ CartesianIndices(E[e_inds...])]
        tlss = [twolevelsystems[pres_inds...][Tuple(i)[2:e_dim]...,pres_inds[e_dim:end]...] for i ∈ CartesianIndices(E[e_inds...])]
        lats = [lattices[pres_inds...][Tuple(i)[2:e_dim]...,pres_inds[e_dim:end]...] for i ∈ CartesianIndices(E[e_inds...])]

        kas = [0 for i ∈ CartesianIndices(E[e_inds...])]
        kbs = [0 for i ∈ CartesianIndices(E[e_inds...])]
        Fs = [[1] for i ∈ CartesianIndices(E[e_inds...])]
        ψ_inits = [[] for i ∈ CartesianIndices(E[e_inds...])]

        if eig_fun == eig_klp
            args = (sims, Es, syss, tlss, lats, kas, kbs, Fs)
        elseif eig_fun == eig_cflp
            ks = [k[e_inds...][i] for i ∈ CartesianIndices(E[e_inds...])]
            args = (sims, ks, syss, tlss, lats, Es, kas, kbs, Fs, ψ_inits)
        else
            throw(ArgumentError("parametric parallization only available for linear eig_k and eig_cf at present."))
        end

        if disp_opt
            T = progress_pmap(eig_fun, args...; progress=Progress(prod(size(Es)),PROGRESS_UPDATE_TIME::Float64,string("parameter dim ", e_dim, " ") ))
        else
            T = pmap(eig_fun, args...)
        end

        for i ∈ CartesianIndices(T)
            E[Tuple(i)[1:e_dim]..., pres_inds[e_dim:end]...] = T[i]
        end
    end
    return E
end

"""
    K = eig_klp(sim, k, systems, twolevelsystems, lattices, ka, kb, F, ψ_init)
"""
function eig_klp(sim::Simulation, k::Number,
    syss::Array{System,M},
    tlss::Array{TwoLevelSystem,M},
    lattices::Array{Bravais,M},
    ka, kb, F, ψ_init) where M

    K = Array{ComplexF64}(undef,length(syss))
    K[1] = k
    for i ∈ 2:length(syss)
        sim = Simulation(sim; sys=syss[i], tls=tlss[i], lattice=lattices[i])
        K[i] = eig_kl(sim, K[i-1], 1, ka, kb, F, ψ_init)[1][1]
    end
    return K
end


"""
    H = eig_cfp(sim, k, systems, twolevelsystems, lattices, η_init, ka, kb, F, ψ_init)
"""
function eig_cfp(sim::Simulation, k::Array{Number,M},
    syss::Array{System,M},
    tlss::Array{TwoLevelSystem,M},
    lattices::Array{Bravais,M},
    η_init, ka, kb, F, ψ_init) where M

    H = Array{ComplexF64}(undef,length(syss))
    H[1] = η_init
    for i ∈ 2:length(syss)
        sim = Simulation(sim; sys=syss[i], tls=tlss[i], lattice=lattices[i])
        H[i] = eig_kl(sim, k[i-1], 1, E[i-1], ka, kb, F, ψ_init)[1][1]
    end
    return K
end

################################################################################
### CONTOUR EIGENSOLVER PARALLEL
################################################################################
"""
    eig_knlp(sim, k, nk, radii; F=[1], r_min=.01, Nq=100, rank_tol=1e-8)
"""
function eig_knlp(sim::Simulation, k::Number, nk::Int, radii, ka=0, kb=0, F=[1], Nq=100, rank_tol=1e-8, r_min=.01)

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
    AA = @distributed (+) for i in 1:Nq

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
        χk′² = sparse(1:N, 1:N, D₀*γ(sim,k′)*F.*F0[:]*k′², N, N)

        A = (∇²+ɛk′²+χk′²)\M
        A₀ = A*dk′/complex(0,2π)
        A₁ = A*k′*dk′/complex(0,2π)

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
            AA₀ = AA*dkk′/complex(0,2π)
            AA₁ = AA*kk′*dkk′/complex(0,2π)

            A₀ = A₀-AA₀
            A₁ = A₁-AA₁
       end
            [A₀ A₁]
    end

    A₀ = AA[:,1:nk]
    A₁ = AA[:,nk .+ (1:nk)]

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
### BAND STRUCTURE
################################################################################

function bands_only_p(
    sim::Simulation,
    k_bloch::Tuple{AbstractArray{S,M},AbstractArray{T,M}},
    num_bands,
    pg = []) where S where T where M

    periodic_boundary_weights!(sim)

    ka_bloch = k_bloch[1]
    kb_bloch = k_bloch[2]

    sims = [sim for i ∈ CartesianIndices(ka_bloch)]
    ks = [0.01 for i ∈ CartesianIndices(ka_bloch)]
    num_bandss = [num_bands for i ∈ CartesianIndices(ka_bloch)]
    Fs = [[1] for i ∈ CartesianIndices(ka_bloch)]
    ψ_inits = [[] for i ∈ CartesianIndices(ka_bloch)]
    kas = [ka_bloch[i] for i ∈ CartesianIndices(ka_bloch)]
    kbs = [kb_bloch[i] for i ∈ CartesianIndices(ka_bloch)]

    if !isempty(pg)
        B = progress_pmap(eig_kl, sims, ks, num_bandss, kas, kbs, Fs, ψ_inits; progress=pg)
    else
        B = pmap(eig_kl, sims, ks, num_bandss, kas, kbs, Fs, ψ_inits)
    end
    bands = Array{ComplexF64, M+1}(undef, num_bands, size(ka_bloch)...)
    for i ∈ CartesianIndices(ka_bloch)
        bands[:,i] = real(B[i][1])
    end

    return bands
end

function band_structure_p(sim::Simulation, num_bloch::Int; num_bands=5, zone=:reduced, interpolation=:cubic, calc_type="band structure ")
    bands, gaps = band_structure(sim::Simulation, num_bloch::Int; num_bands=num_bands, zone=zone, parallel=true, interpolation=interpolation, calc_type=calc_type)
    return bands, gaps
end
function band_structure_p(sim::Simulation, num_bloch::Array{Int,1}; num_bands=5, zone=:reduced, interpolation=nothing, line=nothing)
    bands, gaps, ks = band_structure(sim::Simulation, num_bloch::Array{Int,1}; num_bands=num_bands, zone=zone, parallel=true, interpolation=nothing, line=nothing)
    return bands, gaps, ks
end
function band_structure_p(sim::Simulation, k_bloch::Tuple{AbstractArray{S,M},AbstractArray{T,M}};
    num_bands=5, interpolation=:cubic,
    calc_type = "band structure ",
    pg::Progress=Progress(prod(size(k_bloch[1])), PROGRESS_UPDATE_TIME::Float64, calc_type)) where S where T where M

    bands, gaps, bands_itp = band_structure_p(sim, k_bloch; num_bands=num_bands, parallel=true, interpolation=interpolation, calc_type=calc_type, pg=pg)
    return bands, gaps, bands_itp
end
function band_structure_p(sim::Simulation; waveguide, num_bloch=17, num_bands=5, interpolation=:cubic, zone=:half)
    bands, gaps = band_structure(sim; waveguide=waveguide, num_bloch=num_bloch, num_bands=num_bands, parallel=true, interpolation=interpolation, zone=zone)
    return bands, gaps
end


################################################################################
### WAVEGUIDE DISPERSION
################################################################################
function build_dispersion_p!(sim::Simulation;
    num_wg_bands=round(Int,1.5*max(sim.lat.a/sim.lat.b,sim.lat.b/sim.lat.a)),
    num_bloch=17, num_free_bands=2)

    build_dispersion!(sim; num_wg_bands=num_wg_bands, num_bloch=num_bloch, num_free_bands=num_free_bands, parallel=true)
    return nothing
end

function waveguide_dispersion_p(sim::Simulation, waveguide::Int;
    num_wg_bands=round(Int,1.5*max(sim.lat.a/sim.lat.b,sim.lat.b/sim.lat.a)),
    num_free_bands=2, num_bloch=17, interpolation=:cubic)

    wg_dispersion, bands, gaps, ks = waveguide_dispersion_p(sim, waveguide; num_wg_bands=num_wg_bands, num_free_bands=num_free_bands, num_bloch=num_bloch, parallel=true, interpolation=interpolation)
    return wg_dispersion, bands, gaps, ks
end

################################################################################
### S-MATRIX
################################################################################
"""
    S, flux_sct, flux_tot, absorption = smatrix_p(sim, k; F=[1],
        is_linear=true, disp_opt=true, file_name="",
        num_wg_bands=15, num_bloch=17, num_free_bands=2,
        num_channel_blocks=1)

parallel computation of scattering matrix `S` at frequences `k` on `channels`

`num_channel_blocks` parallelizes the channes into that number of blocks. because
scattering operator can be factorized per frequency and used across many channels,
it is generally not recommended to have this number be anything other than 1.
It is envisioned that this might be useful when only a few frequencies are needed
but many channels (such as might be the case with angular momentum channels).
"""
function smatrix_p(sim::Simulation, k; F=[1],
            is_linear=true, disp_opt=true, file_name="",
            num_wg_bands=round(Int,1.5*max(sim.lat.a/sim.lat.b,sim.lat.b/sim.lat.a)),
            num_bloch=17, num_free_bands=2,
            num_channel_blocks=1)

    build_dispersion_p!(sim; num_wg_bands=num_wg_bands, num_bloch=num_bloch, num_free_bands=num_free_bands)

    # set up iterations
    nc = length(sim.sct.channels)
    nk = length(k)
    channel_inds = vcat(floor.(Int,1:nc/(num_channel_blocks):nc), nc+1)
    num_channel_blocks = length(channel_inds)-1

    sims = [deepcopy(sim) for i ∈ k for j ∈ 1:num_channel_blocks]
    ks = [i for i ∈ k for j ∈ 1:num_channel_blocks]
    channels = [channel_inds[j]:channel_inds[j+1]-1 for i ∈ k for j ∈ 1:num_channel_blocks]
    Fs = [F for i ∈ k for j ∈ 1:num_channel_blocks]
    disp_opts = [false for i ∈ k for j ∈ 1:num_channel_blocks]
    file_names = ["" for i ∈ k for j ∈ 1:num_channel_blocks]

    # parallel map
    if disp_opt
        C = progress_pmap(smatrix_l, sims, ks, channels, Fs, disp_opts, file_names; progress=Progress(length(ks), PROGRESS_UPDATE_TIME::Float64,"smatrix_p "))
    else
        C = pmap(smatrix_l, sims, ks, channels, Fs, disp_opts, file_names)
    end

    # repackage results
    k_inds = [i for i ∈ eachindex(k) for j ∈ 1:num_channel_blocks]
    S = Array{ComplexF64}(undef, nk, nc, nc)
    scattered_flux = Array{Float64}(undef, nk, nc)
    total_flux = Array{Float64}(undef, nk, nc)
    absorption = Array{Float64}(undef, nk, nc)
    for i ∈ eachindex(C)
        S[k_inds[i], channels[i], :] = C[i][1]
        scattered_flux[k_inds[i], channels[i]] = C[i][2]
        total_flux[k_inds[i], channels[i]] = C[i][3]
        absorption[k_inds[i], channels[i]] = C[i][4]
    end

    # save if requested
    if !isempty(file_name)
        fid = open(file_name,"w")
        serialize(fid, (S, flux_sct, flux_tot, absorption, sim, 1))
        close(fid)
    end

    return S, scattered_flux, total_flux, absorption
end
