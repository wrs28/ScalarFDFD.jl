#TODO: Domain: complete documentation for Domain struct
#TODO: Add argument warnings or errors or correcition for some range of alphas and betas so that band scheme can work
"""
    lattice = Bravais(;a=Inf, α=0, b=Inf, β=π/2, x0=0, y0=0)

bravais lattice with:

- `a` is the length of the first vector (by default aligned along x-axis).

- `b` is length of second vector.

- `α` is the angle of first primitive vector.

- `β` is the angle of second primitive vector

- `x0`, `y0` define the origin of the lattice.

"""
struct Bravais
    a::Float64
    b::Float64
    α::Float64
    β::Float64
    x0::Float64
    y0::Float64

    sin²θ::Float64
    cosθ::Float64
    R::Array{Float64,2}
    Rᵀ::Array{Float64,2}
    v1::Array{Float64,1}
    v2::Array{Float64,1}

    function Bravais(;a=Inf, α=0, b=Inf, β=π/2, x0=0, y0=0)
        θ = β-α
        if iszero(θ)
            throw(ArgumentError("lattice angle θ=$(θ) cannot be zero"))
        end
        sinθ, cosθ = sincos(θ)
        R = [ cos(α) -sin(α);
              sin(α)  cos(α)]
        Rᵀ = transpose(R)
        v1 = R*[1,0]
        v2 = R*[cosθ,sinθ]
        new(float(a), float(b), float(α), float(β), float(x0), float(y0),
            sinθ^2, cosθ, R, Rᵀ, v1, v2)
    end
end


"""
    domain = Domain(;is_in_domain=whole_domain, n₁=1, n₂=0, F=0,
        domain_params=[], domain_type=:generic, lattice=Bravais(),
        is_in_region=Function[], region_params=[[]...], n₁_val=[0...],
        n₁_idx=[...], n₂_val=[0...], n₂_idx[...], F_val=[0...], F_idx[...]
        which_asymptote=:none, which_waveguide=0, n₁_val=[1], n₁_idx[1],

where [...] indicates the same size as `n₁_idx`, etc.

- `is_in_domain` is boolean function with inputs `x`, `y`, `domain_params`

- `domain_params` array of floats passed to function is_in_domain

- `n₁`, `n₂` default real and imaginary parts of the index

- `F` the default pump associated with domain

`domain_type`

`lattice` is a bravais lattice

INCOMPLETE DOCUMENTATION

"""
struct Domain
    is_in_domain::Function

    n₁::Float64
    n₂::Float64
    F::Float64
    domain_params::Array{Float64,1}
    domain_type::Symbol
    lattice::Bravais

    is_in_region::Array{Function,1}
    region_params::Array{Array{Float64,1},1}

    n₁_val::Array{Float64,1}
    n₁_idx::Array{Int,1}
    n₂_val::Array{Float64,1}
    n₂_idx::Array{Int,1}
    F_val::Array{Float64,1}
    F_idx::Array{Int,1}

    which_asymptote::Symbol
    which_waveguide::Int

    num_regions::Int

    function Domain(;is_in_domain=whole_domain, n₁=1, n₂=0, F=0, domain_params=[],
        domain_type=:background, lattice=Bravais(),
        is_in_region::Array{Function,1}=Function[],
        region_params=fill(Float64[],length(is_in_region)),
        n₁_val=fill(1,length(is_in_region)), n₁_idx=Array(1:length(n₁_val)),
        n₂_val=fill(0,length(n₁_val)), n₂_idx=Array(1:length(n₂_val)),
        F_val=fill(0,length(n₁_val)), F_idx=Array(1:length(F_val)),
        which_asymptote=:none, which_waveguide=0)

        # argument check
        (
        length(F_idx) == length(n₁_idx) == length(n₂_idx) ==
        length(is_in_region) == length(region_params) ? nothing :
        throw(ArgumentError(
        "lengths of regions don't match.
        F_idx, n₁_idx, n₂_idx, is_in_region, region_params must share length."))
        )

        num_regions= 1 + length(is_in_region)

        if which_asymptote ==:none && which_waveguide !== 0
            throw(ArgumentError(
            "assigned a waveguide identifier $(which_waveguide) to a non-asymptotic domain"))
        end

        return new(is_in_domain, float(n₁), float(n₂), float(F), float.(domain_params),
                    domain_type, lattice, is_in_region, float.(region_params),
                    float.(n₁_val), n₁_idx, float.(n₂_val), n₂_idx, float.(F_val),
                    F_idx, which_asymptote, which_waveguide, num_regions)
    end
end


"""
    sys = System(domains::Array{Domain})
    sys = System(domain)

collection of domains that defines System, an input for Simulation.

input can be array or scalar.
"""
struct System
    domains::Array{Domain,1}

    n_by_region::Array{ComplexF64,1}
    ε_by_region::Array{ComplexF64,1}
    F_by_region::Array{Float64,1}
    domain_by_region::Array{Int,1}

    ε::Array{ComplexF64,2}
    Σ::Array{Array{ComplexF64,1},2}
    F::Array{Float64,2}

    regions::Array{Int,2}
    num_prev_regions::Array{Int,1}

    waveguides::Array{Int,1}

    function System(domains::Array{Domain,1},
                    ε=Array{ComplexF64}(undef, 0, 0),
                    Σ=Array{Array{ComplexF64,1}}(undef, 0, 0),
                    F=Array{Float64}(undef, 0, 0),
                    regions=Array{Int}(undef, 0, 0))

        domains = deepcopy(domains)

        sort!(domains, by=isbulkwaveguide, rev=false)
        sort!(domains, by=isbackground, rev=false)
        sort!(domains, by=isdefect, rev=true)
        sort!(domains, by=ispc, rev=false)
        sort!(domains, by=iswaveguide, rev=true)

        num_regions = 0
        for d ∈ domains
            num_regions += d.num_regions
        end

        n_by_region = Array{ComplexF64}(undef, num_regions)
        F_by_region = Array{Float64}(undef, num_regions)
        domain_by_region = Array{Int}(undef, num_regions)

        num_prev_regions = Array{Int}(undef, length(domains))
        num_prev_regions[1] = 0

        idx = 1
        for i ∈ eachindex(domains)
            d = domains[i]
            if i > 1
                num_prev_regions[i] = domains[i-1].num_regions
            end
            for j ∈ 1:d.num_regions-1
                n_by_region[idx] = complex.(d.n₁_val[d.n₁_idx[j]], d.n₂_val[d.n₂_idx[j]])
                F_by_region[idx] = d.F_val[d.F_idx[j]]
                domain_by_region[idx] = i
                idx += 1
            end
            n_by_region[idx] = complex(d.n₁, d.n₂)
            F_by_region[idx] = d.F
            domain_by_region[idx] = i
            idx += 1
        end
        ε_by_region = n_by_region.^2
        cumsum!(num_prev_regions, num_prev_regions, dims=1)

        waveguides = sort(unique([domains[i].which_waveguide for i ∈ eachindex(domains)]))[2:end]

        return new(domains, n_by_region, ε_by_region, F_by_region,
            domain_by_region, ε, Σ, F, regions, num_prev_regions, waveguides)
    end
end


struct Discretization
    dx::Array{Float64,1}
    sub_pixel_num::Int
    origin::Array{Float64,1}

    N_tr::Array{Int,1}
    N::Array{Int,1}
    dN::Array{Int,2}

    x::Array{Array{Float64,2},1}
    x_tr::Array{Array{Float64,2},1}
    x_idx::Array{Array{Int,1},1}

    X_idx::Array{Int,1}

    function Discretization(dx, sub_pixel_num, origin, N, dN)

        dx = deepcopy(dx)
        sub_pixel_num=deepcopy(sub_pixel_num)
        origin = deepcopy(origin)
        N = deepcopy(N)
        dN = deepcopy(dN)

        N_tr = Array{Int}(undef,2)
        for j ∈ eachindex(N)
            N_tr[j] = N[j] - dN[1,j] - dN[2,j]
        end

        x_tr = [Array{Float64}(undef,N_tr[1],1), Array{Float64}(undef,1,N_tr[2])]
        x = [Array{Float64}(undef,N[1],1), Array{Float64}(undef,1,N[2])]
        x_idx = [Array{Int}(undef, N_tr[1]), Array{Int}(undef, N_tr[2])]
        for j ∈ eachindex(N)
            if !isinf(dx[j])
                x[j][:] .= origin[j] .+ dx[j]*(1/2:N[j]-1/2)
                x_tr[j][:] .= x[j][1] .+ dN[1,j]*dx[j] .+ collect(0:N_tr[j]-1).*dx[j]
                x_idx[j][:] .= dN[1,j] .+ (1:N_tr[j])
            else
                x[j][:] .= origin[j]
                x_tr[j][:] .= x[j][1]
                x_idx[j][:] .= dN[1,j] .+ (1:N_tr[j])
            end
        end

        X_idx = LinearIndices(Array{Bool}(undef, N...))[x_idx...][:]

        return new(float.(dx), sub_pixel_num, float.(origin), N_tr, N, dN,
            x, x_tr, x_idx, X_idx)
    end
end


struct Boundary
    ∂Ω::Array{Float64,2}
    ∂Ω_tr::Array{Float64,2}

    bc::Array{Symbol,2}

    bl::Array{Symbol,2}
    bl_depth::Array{Float64,2}

    indices::Array{Array{Int},1}
    weights::Array{Array{Float64},1}
    shifts::Array{Array{Int},1}

    function Boundary(∂Ω::Array{S,2}, bc::Array{Symbol,2}, bl::Array{Symbol,2},
        bl_depth::Array{T,2}) where S where T

        ∂Ω = deepcopy(∂Ω)
        bc = deepcopy(bc)
        bl = deepcopy(bl)
        bl_depth = deepcopy(bl_depth)

        ∂Ω_tr = Array{Float64}(undef, 2, 2)

        fix_bc!(bc)
        fix_bl!(bl)

        for i ∈ CartesianIndices(bl)
            if bl[i] == :none
                if !iszero(bl_depth[i])
                    @warn "boundary layer [$(i[1]),$(i[2])] set to :none, but has thickness $(bl_depth[i]) > 0
        setting bl_depth[$(i[1]),$(i[2])] = 0"
                end
                bl_depth[i] = 0
            elseif bl[i] ∈ [:pml_out, :pml_in, :abs_out, :abs_in] && bc[i] == :p
                nothing
            elseif bl[i] ∈ [:pml_out, :pml_in, :abs_out, :abs_in] && bc[i] !== :d
                @warn "absorbing layer with non-dirichlet condition for side [$(i[1]), $(i[2])], which is inconsistent
        prioritizing absorbing layer
        setting bc[$(i[1]),$(i[2])] = :d"
                bc[i] = :d
            end
        end

        for j ∈ 1:2
            ∂Ω_tr[1,j] = ∂Ω[1,j] + bl_depth[1,j]
            ∂Ω_tr[2,j] = ∂Ω[2,j] - bl_depth[2,j]
        end

        indices = Array{Array{Int}}(undef, 4)
        weights = Array{Array{Float64}}(undef, 2)
        shifts = Array{Array{Int}}(undef, 4)
        return new(float.(∂Ω), float.(∂Ω_tr), bc, bl, float.(bl_depth), indices, weights, shifts)
    end
end


"""
    channel = Channels(waveguide=0, quantum_number=0)

channel object for Scattering.

`waveguide` identifies which waveguide the channel is defined by

`quantum_number` identifies which mode of the waveguide defines the channel
"""
struct Channels
    waveguide::Int
    quantum_number::Int
    dispersion::Array{AbstractInterpolation,1}
    gaps::Array{Array{Float64,1},1}

    function Channels(waveguide=0, quantum_number=0)
        return  new(waveguide, quantum_number, AbstractInterpolation[], Array{Float64,1}[])
    end
end


"""
    sct = Scattering(channels::Array{Channel})

scattering object for Simulation

`channels` array of Channel objects
"""
struct Scattering
    channels::Array{Channels,1}
    waveguides_used::Array{Int,1}

    ε₀::Array{Array{ComplexF64,2},1}

    function Scattering(channels::Array{Channels,1}=Channels[], ε₀=[])

        channels = deepcopy(channels)
        ε₀ = deepcopy(ε₀)

        num_channels = length(channels)
        waveguides = Array{Int}(undef,num_channels)
        for i ∈ eachindex(waveguides)
            waveguides[i]=channels[i].waveguide
        end
        waveguides_used = unique(waveguides)
        ε₀ = Array{Array{ComplexF64,2},1}(undef,length(waveguides_used))
        return new(channels, waveguides_used, ε₀)
    end
end


"""
    tls = TwoLevelSystem(D₀=0, k₀=10, γp=1e8, ω₀=k₀)

`D₀` is pump parameter

`k₀` is atomic frequency

`γp` is depolarization rate
"""
struct TwoLevelSystem
    D₀::Float64
    k₀::Float64
    γp::Float64
    ω₀::Float64

    function TwoLevelSystem(D₀=0., k₀=10., γp=1e8, ω₀=k₀)
            ω₀ = k₀
        return new(float(D₀), float(k₀), float(γp), float(ω₀))
    end
end


"""
See also: [`System`](@ref), [`Boundary`](@ref), [`Discretization`](@ref), [`Scattering`](@ref), [`TwoLevelSystem`](@ref), [`Bravais`](@ref)

    sim = Simulation(; sys, bnd, dis, sct=Scattering(), tls=TwoLevelSystem(), lat=Bravais(bnd))

simulation object
"""
struct Simulation
    sys::System
    bnd::Boundary
    dis::Discretization
    sct::Scattering
    tls::TwoLevelSystem
    lat::Bravais

    function Simulation(;sys::System=System(Domain()), bnd::Boundary,
        dis::Discretization,
        sct::Scattering=Scattering(),
        tls::TwoLevelSystem=TwoLevelSystem(),
        lat::Bravais=Bravais(bnd),
        disp_opt=false)

        sys = deepcopy(sys)
        bnd = deepcopy(bnd)
        dis = deepcopy(dis)
        sct = deepcopy(sct)
        tls = deepcopy(tls)
        lat = deepcopy(lat)

        ∂Ω = bnd.∂Ω
        L = Array{Float64}(undef,2)
        for i ∈ 1:2
            L[i] = ∂Ω[2,i]-∂Ω[1,i]
        end

        dx = dis.dx
        try
            N = round.(Int,L./dx)
        catch
            throw(ArgumentError("lattice spacings $dx inconsistent with size $L"))
        end

        dN = Array{Int}(undef,2,2)
        index1 = findmin(N)[2]
        index2 = mod1(index1+1,2)
        if N[index1]≤1
            N[index1] = 1
            dx[index1] = Inf
            dx[index2] = L[index2]/N[index2]
            dN[:,index1] .= 0
            dN[:,index2] .= floor.(Int, bnd.bl_depth[:,index2]/dx[index2])
            bnd.bl_depth[:,index1] .= 0
            bnd.bl_depth[:,index2] .= dN[:,index2]*dx[index2]
            bnd.bl[:,index1] .= :none
            bnd.bc[:,index1] .= :d
        else
            dx[index1] = L[index1]/N[index1]
            dx[index2] = L[index2]/round(Int,L[index2]/dx[index1])
            N = round.(Int,L./dx)
            for i ∈ 1:2
                dN[:,i] = floor.(Int, bnd.bl_depth[:,i]/dx[i])
                bnd.bl_depth[:,i] = dN[:,i]*dx[i]
            end
        end

        origin = bnd.∂Ω[1,:]

        bnd = Boundary(∂Ω, bnd.bc, bnd.bl, bnd.bl_depth)
        dis = Discretization(dx, dis.sub_pixel_num, origin, N, dN)

        if :p ∈ bnd.bc[:,1] && isinf(lat.a)
            @warn "periodic bc in dimension 1, but infinite lattice constant. redefining lat.a = width of domain"
            lat.a = bnd.∂Ω[2,1]-bnd.∂Ω[1,1]
        end
        if :p ∈ bnd.bc[:,2] && isinf(lat.b)
            @warn "periodic bc in dimension 2, but infinite lattice constant. redefining lat.b = height of domain"
            lat.b = bnd.∂Ω[2,2]-bnd.∂Ω[1,2]
        end
        for i ∈ eachindex(sys.domains)
            if N[1] == 1 && !isinf(sys.domains[i].lattice.a)
                @warn "one-dimensional (vertical) system, but domain $(i) has finite transverse lattice constant $(sys.domains[i].lattice.a) < ∞"
            end
            if N[2] == 1 && !isinf(sys.domains[i].lattice.b)
                @warn "one-dimensional (horizontal) system, but domain $(i) has finite transverse lattice constant $(sys.domains[i].lattice.b) < ∞"
            end
        end

        ɛ, F, regions = ScalarFDFD.sub_pixel_smoothing(bnd, dis, sys; disp_opt=disp_opt)
        Σ = ScalarFDFD.σ(bnd, dis)
        sys = System(sys.domains, ε, Σ, F, regions)

        for i ∈ 1:length(sct.waveguides_used)
            sct.ɛ₀[i] = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
            wg_inds = findall([sys.domains[j].which_waveguide for j ∈ eachindex(sys.domains)] .== sct.waveguides_used[i])
            wg_domains = Array{Domain}(undef,length(wg_inds))
            for j ∈ eachindex(wg_domains)
                wg_domains[j] = Domain(sys.domains[wg_inds[j]]; :which_asymptote => :none, :which_waveguide => 0)
            end
            if isempty(wg_domains)
                @warn "channels references undefined waveguide $(sct.waveguides_used[i]), scattering calculations will fail"
            else
                sct.ɛ₀[i][:] = ScalarFDFD.sub_pixel_smoothing(bnd, dis, System(wg_domains), disp_opt=disp_opt)[1][:]
            end
        end

        return new(sys, bnd, dis, sct, tls, lat)
    end
end
