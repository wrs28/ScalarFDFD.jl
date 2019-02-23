################################################################################
########## BAND STRUCTURE
################################################################################
struct Mnemosyne end

"""
    p = plot(sim, bands, gaps, kxy; options...)
    p = plot(sim, bands, kxy; options...)
    p = plot(sim, bands, gaps; options...)
    p = plot(sim, bands; options...)

plot band structure as computed from `band_structure`

options list:

`universal_scale`=1, which divides all frequencies.

`band_width`=2, `band_style`=:solid

`gap_color`=:lightgrey, `gap_width`=1, `gap_style`=:dash

`n`=1:size(bands,1) which bands to plot. If scalar, plot is false color plot.
If vector, then surface plots.

note everything between `sim` and `ks` can appear in any order.
"""
@recipe function f(::Mnemosyne, sim::Simulation,
    bands::Array{Float64,N},
    gaps::Tuple{Array{Float64,1},Array{Float64,1}},
    ks::AbstractArray,
    n,
    universal_scale::Number,
    band_width, band_style::Symbol,
    gap_color, gap_width, gap_style) where N

    title --> "Band Structure"
    legend --> false
    grid --> false

    ks = ks*universal_scale
    bands = bands*universal_scale
    gaps = (gaps[1]*universal_scale, gaps[2]*universal_scale)

    if N == 2
        if universal_scale==1
            xlabel --> LaTeXString("propagation constant  \\beta")
            ylabel --> LaTeXString("frequency  \\omega")
        else
            xlabel --> string("propagation constant ", LaTeXString("\\beta"), "/2", LaTeXString("\\pi"), fmt("2.2f", 2π*universal_scale))
            ylabel --> string("frequency ", LaTeXString("\\omega"), "/2", LaTeXString("\\pi"), fmt("2.2f", 2π*universal_scale))
        end

        if !isempty(ks)
            @series begin
                linewidth --> band_width
                linestyle --> band_style
                ks, bands'
            end
            @series begin
                linecolor := gap_color
                linestyle --> gap_style
                linewidth --> gap_width
                [ks[1];ks[end]], [gaps[1] gaps[1]]'
            end
            @series begin
                linecolor := gap_color
                linestyle --> gap_style
                linewidth --> gap_width
                [ks[1];ks[end]], [gaps[2] gaps[2]]'
            end
        else
            xlabel := ""
            xticks := []
            @series begin
                linewidth --> band_width
                linestyle --> band_style
                bands'
            end
            @series begin
                linecolor := gap_color
                linestyle --> gap_style
                linewidth --> gap_width
                [1;size(bands,2)], [gaps[1] gaps[1]]'
            end
            @series begin
                linecolor := gap_color
                linestyle --> gap_style
                linewidth --> gap_width
                [1;size(bands,2)], [gaps[2] gaps[2]]'
            end
        end
    elseif N > 2
        if !(typeof(n) <: Union{AbstractArray{Int},Int})
            throw(ArgumentError("argument `n`=$n is neither an array nor an int"))
        end
        ylabel --> LaTeXString("\\omega")
        if !isempty(size(n))
            if !isempty(ks)
                xlabel --> LaTeXString("k")
                for i ∈ n
                    @series begin
                        seriestype := :surface
                        ks[1],ks[2], bands[i,:,:]'
                    end
                end
            else
                for i ∈ n
                    @series begin
                        seriestype := :surface
                        bands[i,:,:]'
                    end
                end
            end
        else
            if !isempty(ks)
                xlabel --> LaTeXString("k")
                @series begin
                    seriestype --> :contour
                    fill --> true
                    colorbar --> true
                    ks[1],ks[2],bands[n,:,:]'
                end
            else
                @series begin
                    seriestype --> :contour
                    fill --> true
                    colorbar --> true
                    bands[n,:,:]'
                end
            end
        end
    else
        nothing
    end
end


@recipe function f(::Mnemosyne,
    sim::Simulation,
    bands::Array{Array{Array{Array{Float64,2},1},1},1},
    gaps::Tuple,
    ks,
    n,
    universal_scale::Number,
    band_width, band_style::Symbol,
    gap_color, gap_width, gap_style)

    num_bands = size(bands[1][1][1],1)
    layout := length(bands)
    grid --> false
    for i ∈ 1:length(bands)
        for j ∈ 1:length(bands[i])
            B = bands[i][j][1]
            zone = Array{Int}(undef,2)
            for m ∈ 2:3
                zone[m-1] = size(B,2)
                B = hcat(B[:,:], bands[i][j][m][:,2:end])
            end
            for m ∈ 1:2
                @series begin
                    legend := false
                    subplot :=i
                    linecolor := gap_color
                    lw := 1.5
                    [zone[m], zone[m]], [0, maximum(B)]
                end
            end
            @series begin
                subplot := i
                (Mnemosyne(), sim, B, gaps, Float64[], n, universal_scale, band_width, band_style, gap_color, gap_width, gap_style)
            end
        end
    end
end

@recipe function f(sim::Simulation,
    bands::Union{Array{Float64,N}, Array{Array{Array{Array{Float64,2},1},1},1}},
    args...;
    n = 1:size(bands,1),
    universal_scale = 1,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE) where N

    outargs = outargs_for_band_plotting(args)
    (Mnemosyne(), sim, bands, outargs..., n, universal_scale, band_width, band_style, gap_color, gap_width, gap_style)
end
@recipe function f(sim::Simulation,
    arg1,
    bands::Union{Array{Float64,N},Array{Array{Array{Array{Float64,2},1},1},1}},
    args...;
    n = 1:size(bands,1),
    universal_scale = 1,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE) where N

    outargs = outargs_for_band_plotting((arg1,args...))
    (Mnemosyne(), sim, bands, outargs..., n, universal_scale, band_width, band_style, gap_color, gap_width, gap_style)
end
@recipe function f(sim::Simulation,
    arg1,
    arg2,
    bands::Union{Array{Float64,N},Array{Array{Array{Array{Float64,2},1},1},1}};
    n = 1:size(bands,1),
    universal_scale = 1,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE) where N

    outargs = outargs_for_band_plotting((arg1,arg2))
    (Mnemosyne(), sim, bands, outargs..., n, universal_scale, band_width, band_style, gap_color, gap_width, gap_style)
end


"""
    outargs = outargs_for_band_plotting
"""
function outargs_for_band_plotting(args)

    outargs = Array{Any}(undef,2)
    outargs[1] = (Float64[],Float64[])
    outargs[2] = []
    for a ∈ args
        if typeof(a) <: Tuple
            outargs[1] = a
        elseif typeof(a) <: AbstractArray
            outargs[2] = a
        end
    end

    return outargs
end



################################################################################
########## WG DISPERSION
################################################################################


struct Iapetos end

"""
    p = plot(sim, bands, gaps, wg_dispersion, ks; options...)
    p = plot(sim, bands, wg_dispersion, ks; options...)
    p = plot(sim, gaps, wg_dispersion, ks; options...)
    p = plot(sim, ks, wg_dispersion; options...)
    p = plot(sim, wg_dispersion; options...)

options list:

`band_color`=:darkgrey, `band_width`=1.2, `band_style`=:solid

`gap_color`=:lightgrey, `gap_width`=1, `gap_style`=:dash

`dispersion_width`=3, `dispersion_style`=:solid

note everything after `sim` can appear in any order
"""
@recipe function f(::Iapetos,
    sim::Simulation,
    wg_dispersion::Array{Interpolations.AbstractInterpolation},
    bands::Array{Float64,2},
    gaps::Tuple,
    ks::AbstractArray{Float64,1},
    band_color,
    band_width,
    band_style,
    gap_color,
    gap_width,
    gap_style,
    dispersion_width,
    dispersion_style)

    legend --> false
    title --> "Band Structure and Waveguide Dispersion"
    xlabel --> LaTeXString("propagation constant  \\beta")
    ylabel --> LaTeXString("frequency  \\omega")
    grid --> false

    if isempty(ks)
        ks = wg_dispersion[1].ranges[1]
    end
    for i ∈ 1:length(wg_dispersion)
        kwg = wg_dispersion[i].ranges[1][1]:wg_dispersion[i].ranges[1][end]/201:wg_dispersion[i].ranges[1][end]
        @series begin
            linewidth --> dispersion_width
            linestyle --> dispersion_style
            kwg, wg_dispersion[i].(kwg)
        end
    end
    @series begin
        color := band_color
        linewidth --> band_width
        linestyle --> band_style
        ks, bands'
    end
    @series begin
        color := gap_color
        linewidth --> gap_width
        linestyle --> gap_style
        [ks[1];ks[end]], [gaps[1] gaps[1]]'
    end
    @series begin
        color := gap_color
        linewidth --> gap_width
        linestyle --> gap_style
        [ks[1];ks[end]], [gaps[2] gaps[2]]'
    end
end

@recipe function f(sim::Simulation,
    wg_dispersion::Array{Interpolations.AbstractInterpolation,1},
    args...;
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    outargs = outargs_for_wg_plotting(args)

    (Iapetos(), sim, wg_dispersion, outargs..., band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end
@recipe function f(sim::Simulation,
    arg1,
    wg_dispersion::Array{Interpolations.AbstractInterpolation,1},
    args...;
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    args = [args..., arg1]
    outargs = outargs_for_wg_plotting(args)

    (Iapetos(), sim, wg_dispersion, outargs..., band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end
@recipe function f(sim::Simulation,
    arg1, arg2,
    wg_dispersion::Array{Interpolations.AbstractInterpolation,1},
    args...;
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    args = [args..., arg1, arg2]
    outargs = outargs_for_wg_plotting(args)

    (Iapetos(), sim, wg_dispersion, outargs..., band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end
@recipe function f(sim::Simulation,
    arg1, arg2, arg3,
    wg_dispersion::Array{Interpolations.AbstractInterpolation,1};
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    args = [arg1, arg2, arg3]
    outargs = outargs_for_wg_plotting(args)

    (Iapetos(), sim, wg_dispersion, outargs..., band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end
@recipe function f(simd::Simulation, sim::Simulation;
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    wg_dispersion = AbstractInterpolation[sim.sct.channels[i].dispersion[1] for i ∈ eachindex(sim.sct.channels)]
    gaps = ([sim.sct.channels[i].gaps[j][1] for i ∈ eachindex(sim.sct.channels) for j ∈ eachindex(sim.sct.channels[i].gaps)],[sim.sct.channels[i].gaps[j][2] for i ∈ eachindex(sim.sct.channels) for j ∈ eachindex(sim.sct.channels[i].gaps)])

    (Iapetos(), sim, wg_dispersion, Array{Float64}(undef,0,0), gaps, Float64[], band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end
@recipe function f(sim::Simulation, sct::Scattering;
    band_color = BAND_COLOR,
    band_width = BAND_WIDTH,
    band_style = BAND_STYLE,
    gap_color = GAP_COLOR,
    gap_width = GAP_WIDTH,
    gap_style = GAP_STYLE,
    dispersion_width = DISPERSION_WIDTH,
    dispersion_style = DISPERSION_STYLE)

    wg_dispersion = AbstractInterpolation[sct.channels[i].dispersion[1] for i ∈ eachindex(sct.channels)]
    gaps = ([sct.channels[i].gaps[j][1] for i ∈ eachindex(sim.sct.channels) for j ∈ eachindex(sct.channels[i].gaps)],[sct.channels[i].gaps[j][2] for i ∈ eachindex(sct.channels) for j ∈ eachindex(sct.channels[i].gaps)])

    (Iapetos(), sim, wg_dispersion, Array{Float64}(undef,0,0), gaps, Float64[], band_color, band_width, band_style, gap_color, gap_width, gap_style, dispersion_width, dispersion_style)
end

"""
    outargs = outargs_for_wg_plotting(args)
"""
function outargs_for_wg_plotting(args)
    outargs = Array{Any}(undef,3)
    outargs[1] = Array{Float64}(undef,0,0)
    outargs[2] = (Float64[],Float64[])
    outargs[3] = Float64[]
    for a ∈ args
        if typeof(a) <: Array{Float64,2}
            outargs[1] = a
        elseif typeof(a) <: Tuple
            outargs[2] = a
        elseif typeof(a) <: AbstractArray{Float64,1}
            outargs[3] = a
        end
    end
    return outargs
end
