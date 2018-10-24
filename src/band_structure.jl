#TODO : animate bands

################################################################################
### WRAPPERS
################################################################################
"""
    bands, gaps = band_structure(sim, num_bloch::Int; num_bands=5, zone=:half, interpolation=:cubic)

band structure along path specified by `zone`, which must be one of :reduced, :half, :full,
or an integer or array of integers.

Zone specification is

| zone | number |
|:----:|:----:|
| 3    |  4   |
| 1    |  2   |

Except for reduced, each zone is defined by (_ ↑ ←, ↗ ↓ ↖)
"""
function band_structure(sim::Simulation, num_bloch::Int; num_bands=5, zone=:half,
        parallel=nprocs()>1, interpolation=:cubic, calc_type="band structure ",
        disp_opt=true)

    k_paths, zones, order = band_paths(sim, num_bloch, zone)
    bands = Array{Array{Array{Array{Float64,2},1},1},1}(undef, length(zones))
    if disp_opt
        pg = Progress(length(zones)*length(zones[1])*3*num_bloch, PROGRESS_UPDATE_TIME::Float64, calc_type)
    else
        pg = nothing
    end
    for i ∈ eachindex(bands)
        bands[i] = Array{Array{Array{Float64,2},1},1}(undef, length(zones[i]))
        for j ∈ eachindex(zones[i])
            bands[i][j] = Array{Array{Float64,2},1}(undef, 3)
            for m ∈ 1:3
                kb = (k_paths[zones[i][j][m]][1][order[i][j][m]], k_paths[zones[i][j][m]][2][order[i][j][m]])
                if parallel
                    bands[i][j][m] = real(bands_only_p(sim, kb, num_bands, pg))
                else
                    bands[i][j][m] = real(bands_only(sim, kb, num_bands, pg))
                end
            end
        end
    end
    gaps, bands_itp = find_gaps(bands, interpolation)

    return bands, gaps
end


"""
    bands, gaps, kxy = band_structure(sim, num_bloch::Array; num_bands=5, zone=:half)

band structure as full surface over zone specified by `zone`, which must be one of
:reduced, :half, :full.
"""
function band_structure(sim::Simulation, num_bloch::Array{Int,1}; num_bands=5, zone=:half,
    parallel=nprocs()>1, calc_type="band structure ", disp_opt=true, interpolation=nothing,
    line=nothing)

    if length(num_bloch)==1
        num_bloch = fill(num_bloch[1],2)
    end

    if zone == :reduced
        kx = (0:1/num_bloch[1]:1)*2π/2sim.lat.a
        ky = (0:1/num_bloch[2]:1)*2π*sqrt(sim.lat.sin²θ)/2sim.lat.b
    elseif zone == :half
        kx = (0:1/num_bloch[1]:1)*2π/sim.lat.a
        ky = (0:1/num_bloch[2]:1)*2π*sqrt(sim.lat.sin²θ)/2sim.lat.b
    elseif zone == :full
        kx = (0:1/num_bloch[1]:1)*2π/sim.lat.a
        ky = (0:1/num_bloch[2]:1)*2π*sqrt(sim.lat.sin²θ)/sim.lat.b
    else
        throw(ArgumentError("invalid zone $zone, must be one of :reduced, :half, :full"))
    end
    KAB = bravais_coordinates.(kx,ky',Ref(Bravais(sim.lat; :a=>2π/sim.lat.a, :b=>2π/sim.lat.b)))
    ka = Array{Float64}(undef, size(KAB))
    kb = Array{Float64}(undef, size(KAB))
    for i ∈ eachindex(KAB)
        ka[i] = KAB[i][1]
        kb[i] = KAB[i][2]
    end

    if disp_opt
        pg = Progress(prod(size(ka)), PROGRESS_UPDATE_TIME, calc_type)
    else
        pg = nothing
    end
    if parallel
        bands = bands_only_p(sim, (ka, kb), num_bands, pg)
    else
        bands = bands_only(sim, (ka, kb), num_bands, pg)
    end
    gaps, bands_itp = find_gaps([[[real(bands)]]], interpolation)

    return real(bands), gaps, [kx, ky]
end


function band_structure(sim::Simulation, k_bloch::Tuple{AbstractArray{S,M},AbstractArray{T,M}};
    num_bands=5, parallel=nprocs()>1, interpolation=:cubic,
    calc_type = "band structure ",
    pg::Progress=Progress(prod(size(k_bloch[1])), PROGRESS_UPDATE_TIME, calc_type),
    disp_opt=true) where S where T where M

    if !disp_opt
        pg = nothing
    end
    if parallel
        bands = real(bands_only_p(sim, k_bloch, num_bands, pg))
    else
        bands = real(bands_only(sim, k_bloch, num_bands, pg))
    end
    gaps, bands_itp = find_gaps([[[bands]]], interpolation)

    return bands, gaps, bands_itp
end


function band_structure(sim::Simulation, k_bloch::Tuple{AbstractArray{S,M},AbstractArray{T,M}},
    num_bloch_interpolations::Int, reciprocal_basis=true;
    num_bands=5, parallel=nprocs()>1, interpolation=:cubic,
    calc_type = "band structure ",
    pg::Progress=Progress(num_bloch_interpolations*(length(k_bloch[1])-1), PROGRESS_UPDATE_TIME, calc_type),
    disp_opt=true) where S where T where M

    kas = Array{Float64}(undef, num_bloch_interpolations*(length(k_bloch[1])-1))
    kbs = deepcopy(kas)
    for i ∈ 1:length(k_bloch[1])-1
        if reciprocal_basis
            ka_start = k_bloch[1][i]*2π/sim.lat.a
            kb_start = k_bloch[2][i]*2π/sim.lat.b

            ka_stop = k_bloch[1][i+1]*2π/sim.lat.a
            kb_stop = k_bloch[2][i+1]*2π/sim.lat.b
        else
            ka_start = k_bloch[1][i]
            kb_start = k_bloch[2][i]

            ka_stop = k_bloch[1][i+1]
            kb_stop = k_bloch[2][i+1]
        end
        ka_temp = LinRange(ka_start, ka_stop, num_bloch_interpolations)
        kb_temp = LinRange(kb_start, kb_stop, num_bloch_interpolations)
        for j ∈ 1:num_bloch_interpolations
            kas[(i-1)*num_bloch_interpolations + j] = ka_temp[j]
            kbs[(i-1)*num_bloch_interpolations + j] = kb_temp[j]
        end
    end

    bands, gaps, bands_itp = band_structure(sim, (kas, kbs); num_bands=num_bands, parallel=parallel,
        interpolation=interpolation, calc_type=calc_type, pg=pg, disp_opt=disp_opt)

    return bands, gaps, bands_itp, (kas, kbs)
end


"""
    bands, gaps, ks = band_structure(sim; waveguide, num_bloch=17, num_bands=5, interpolation=:cubic)

bands and gaps of photonic crystal that defines background of `waveguide`
"""
function band_structure(sim::Simulation; waveguide, num_bloch=17, num_bands=5,
    parallel=nprocs()>1, interpolation=:cubic, zone=:half, disp_opt=true)

    wg_sim = extract_waveguide_simulation(sim, waveguide)
    pc_sim = Simulation(2, wg_sim)
    bands, gaps = band_structure(pc_sim, num_bloch; num_bands=num_bands,
        parallel=parallel, interpolation=interpolation, zone=zone,
        calc_type="waveguide band structure ", disp_opt=disp_opt)
    return bands, gaps
end

################################################################################
### CORE
################################################################################
"""
    bands = bands_only(sim, k_bloch::Tuple, num_bands, progress)
"""
function bands_only(
    sim::Simulation,
    k_bloch::Tuple{AbstractArray{S,M},AbstractArray{T,M}},
    num_bands,
    pg=nothing) where S where T where M

    ka_bloch = k_bloch[1]
    kb_bloch = k_bloch[2]

    bands = Array{ComplexF64, M+1}(undef, num_bands, size(ka_bloch)...)
    for i ∈ CartesianIndices(ka_bloch)
        bands[:,i] = eig_k(sim, 0.001, num_bands; ka=ka_bloch[i], kb=kb_bloch[i])[1]

        if typeof(pg)==Progress
            next!(pg)
        end
    end

    return bands
end


"""
    gaps, bands_itp = find_gaps(bands, inerpolation)
"""
function find_gaps(bands::Array{Array{Array{Array{Float64,M},1},1},1}, interpolation) where M

    num_bands = size(bands[1][1][1],1)
    bands_itp = Array{AbstractExtrapolation}(undef, num_bands)

    if M==2 # not surface
        band_tops = Array{Float64}(undef, num_bands)
        band_bottoms = Array{Float64}(undef, num_bands)
        band_tops = zeros(Float64, num_bands)
        band_bottoms = Inf*ones(Float64, num_bands)
        for i ∈ eachindex(bands)
            for j ∈ eachindex(bands[i])
                for m ∈ eachindex(bands[i][j])
                    tops = findmax(real(bands[i][j][m]); dims=2)[2]
                    bottoms = findmin(real(bands[i][j][m]); dims=2)[2]
                    for n ∈ 1:num_bands
                        if interpolation == :linear
                            interpolator = BSpline(Linear())
                        elseif interpolation == :quadratic
                            interpolator = BSpline(Quadratic(Line(OnGrid())))
                        elseif interpolation == :cubic
                            interpolator = BSpline(Cubic(Line(OnGrid())))
                        else
                            throw(ArgumentError("invalid order $(interpolation)"))
                        end
                        bands_itp[n] = extrapolate(interpolate(bands[i][j][m][n,:], interpolator), Reflect())
                        band_tops[n] = max(band_tops[n], bands_itp[n](optimize(x->-bands_itp[n](x[1]), [float(tops[n][2])], BFGS()).minimizer)[1])::Float64
                        band_bottoms[n] = min(band_bottoms[n], bands_itp[n](optimize(x->abs(bands_itp[n](x[1])), [float(bottoms[n][2])], BFGS()).minimizer)[1])::Float64
                    end
                end
            end
        end
        gaps = band_bottoms[2:end]-band_tops[1:end-1]
        gaps_start = band_tops[1:end-1][gaps.>0]
        gaps_stop = band_tops[1:end-1][gaps.>0] + gaps[gaps.>0]
    else #surface
        bands = bands[1][1][1]
        tops = (findmax(bands; dims=(2,3))[1])
        bottoms = findmin(bands; dims=(2,3))[1]
        gaps = bottoms[2:end]-tops[1:end-1]
        gaps_start = tops[1:end-1][gaps.>0]
        gaps_stop = tops[1:end-1][gaps.>0] + gaps[gaps.>0]
    end

    return (gaps_start, gaps_stop), bands_itp
end

################################################################################
### AUXILLIARIES
################################################################################
"""
    k_paths, zones, order = band_paths(sim, num_bloch, zone)

`zone` can be an integer, array, or one of `:reduced`, `:half`, `:full`
"""
function band_paths(sim::Simulation, num_bloch, zone=:half)

    k1 = [0,0]
    k2 = 2π*(sim.lat.v1/2sim.lat.a)
    k3 = 2π*(sim.lat.v1/1sim.lat.a)

    k4 = 2π*(sim.lat.v2/2sim.lat.b)
    k5 = 2π*(sim.lat.v1/2sim.lat.a+sim.lat.v2/2sim.lat.b)
    k6 = 2π*(sim.lat.v1/1sim.lat.a+sim.lat.v2/2sim.lat.b)

    k7 = 2π*(sim.lat.v2/1sim.lat.b)
    k8 = 2π*(sim.lat.v1/2sim.lat.a+sim.lat.v2/1sim.lat.b)
    k9 = 2π*(sim.lat.v1/1sim.lat.a+sim.lat.v2/1sim.lat.b)


    ka_bloch01 = LinRange(k1[1], k2[1], num_bloch)
    ka_bloch02 = LinRange(k2[1], k3[1], num_bloch)

    ka_bloch03 = LinRange(k1[1], k4[1], num_bloch)
    ka_bloch04 = LinRange(k2[1], k5[1], num_bloch)

    ka_bloch05 = LinRange(k1[1], k5[1], num_bloch)
    ka_bloch06 = LinRange(k2[1], k4[1], num_bloch)
    ka_bloch07 = LinRange(k2[1], k6[1], num_bloch)
    ka_bloch08 = LinRange(k3[1], k5[1], num_bloch)

    ka_bloch09 = LinRange(k4[1], k5[1], num_bloch)
    ka_bloch10 = LinRange(k5[1], k6[1], num_bloch)

    ka_bloch11 = LinRange(k4[1], k7[1], num_bloch)
    ka_bloch12 = LinRange(k5[1], k8[1], num_bloch)

    ka_bloch13 = LinRange(k4[1], k8[1], num_bloch)
    ka_bloch14 = LinRange(k5[1], k7[1], num_bloch)
    ka_bloch15 = LinRange(k5[1], k9[1], num_bloch)
    ka_bloch16 = LinRange(k6[1], k8[1], num_bloch)


    kb_bloch01 = LinRange(k1[2], k2[2], num_bloch)
    kb_bloch02 = LinRange(k2[2], k3[2], num_bloch)

    kb_bloch03 = LinRange(k1[2], k4[2], num_bloch)
    kb_bloch04 = LinRange(k2[2], k5[2], num_bloch)

    kb_bloch05 = LinRange(k1[2], k5[2], num_bloch)
    kb_bloch06 = LinRange(k2[2], k4[2], num_bloch)
    kb_bloch07 = LinRange(k2[2], k6[2], num_bloch)
    kb_bloch08 = LinRange(k3[2], k5[2], num_bloch)

    kb_bloch09 = LinRange(k4[2], k5[2], num_bloch)
    kb_bloch10 = LinRange(k5[2], k6[2], num_bloch)

    kb_bloch11 = LinRange(k4[2], k7[2], num_bloch)
    kb_bloch12 = LinRange(k5[2], k8[2], num_bloch)

    kb_bloch13 = LinRange(k4[2], k8[2], num_bloch)
    kb_bloch14 = LinRange(k5[2], k7[2], num_bloch)
    kb_bloch15 = LinRange(k5[2], k9[2], num_bloch)
    kb_bloch16 = LinRange(k6[2], k8[2], num_bloch)


    k_bloch01 = (ka_bloch01, kb_bloch01)
    k_bloch02 = (ka_bloch02, kb_bloch02)
    k_bloch03 = (ka_bloch03, kb_bloch03)
    k_bloch04 = (ka_bloch04, kb_bloch04)
    k_bloch05 = (ka_bloch05, kb_bloch05)
    k_bloch06 = (ka_bloch06, kb_bloch06)
    k_bloch07 = (ka_bloch07, kb_bloch07)
    k_bloch08 = (ka_bloch08, kb_bloch08)
    k_bloch09 = (ka_bloch09, kb_bloch09)
    k_bloch10 = (ka_bloch10, kb_bloch10)
    k_bloch11 = (ka_bloch11, kb_bloch11)
    k_bloch12 = (ka_bloch12, kb_bloch12)
    k_bloch13 = (ka_bloch13, kb_bloch13)
    k_bloch14 = (ka_bloch14, kb_bloch14)
    k_bloch15 = (ka_bloch15, kb_bloch15)
    k_bloch16 = (ka_bloch16, kb_bloch16)

    k_paths = (k_bloch01, k_bloch02, k_bloch03, k_bloch04,
        k_bloch05, k_bloch06, k_bloch07, k_bloch08,
        k_bloch09, k_bloch10, k_bloch11, k_bloch12,
        k_bloch13, k_bloch14, k_bloch15, k_bloch16)

    zones = (((1,4,9),(5,3,6)),
             ((2,4,10),(8,3,7)),
             ((1,12,9),(14,11,13)),
             ((2,12,10),(15,11,16)))

    polarity = (((1,1,-1),(-1,1,-1)),
                ((-1,1,1),(-1,1,-1)),
                ((1,-1,-1),(1,-1,1)),
                ((-1,-1,1),(1,-1,1)))
    if zone==:reduced
        zones = (((1,4,5),),)
        polarity = (((1,1,-1),),)
    elseif zone==:half
        zones = (((1,4,9),(5,3,6)),
                 ((2,4,10),(8,3,7)))
        polarity = (((1,1,-1),(-1,1,-1)),
                    ((-1,1,1),(-1,1,-1)))
    elseif zone==:full
        nothing
    elseif typeof(zone) <: Int
        zones = zones[[zone]]
        polarity = polarity[[zone]]
    elseif typeof(zone) <: Array{Int}
        zones = zones[zone]
        polarity = polarity[zone]
    else
        throw(ArgumentError("invalid zone specification $zone, must be one of :reduced, :half, :full, or an integer or array of integers."))
    end

    order = Array{Array{Array{AbstractArray}}}(undef,length(zones))
    for i ∈ 1:length(order)
        order[i] = Array{Array{AbstractArray}}(undef, length(zones[i]))
        for j ∈ 1:length(zones[i])
            order[i][j] = Array{AbstractArray}(undef, 3)
            for m ∈ 1:3
                if polarity[i][j][m]>0
                    step = 1
                    start = 1
                    stop = num_bloch
                else
                    step = -1
                    start = num_bloch
                    stop = 1
                end
                order[i][j][m] = start:step:stop
            end
        end
    end

    return k_paths, zones, order
end
