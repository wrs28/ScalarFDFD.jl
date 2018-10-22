#TODO: check 1d plots, to waveguixe_dispersion

BAND_COLOR = :lightgrey
BAND_WIDTH = 1.2
BAND_STYLE = :solid
GAP_COLOR = :lightgrey
GAP_WIDTH = 1
GAP_STYLE = :dash
DISPERSION_WIDTH = 3
DISPERSION_STYLE = :solid

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

################################################################################
########## SIMULATION
################################################################################
"""
    p = plot(sim, theme=:DEFAULT_COLOR_SCHEME)
"""
@recipe function f(sim::Simulation, theme::Symbol=DEFAULT_COLOR_SCHEME::Symbol)
    (sim, ComplexF64[], theme)
end


################################################################################
########## SOLUTIONS
################################################################################
"""
    p = plot(sim, ψ, theme=:DEFAULT_COLOR_SCHEME; by=nothing, truncate=true, seriestype=:heatmap)

plots
to turn off translucent effect, add optional argument `seriesalpha=0`
vary type of plot with `seriestype`, e.g. `seriestype=:contour`
"""
@recipe function f(sim::Simulation, ψ::Array{ComplexF64}, theme::Symbol=DEFAULT_COLOR_SCHEME::Symbol; by=nothing, truncate=true, seriestype=:heatmap)

    # check whether ψ was originally provided or if just plotting sim
    if length(ψ) > 0
        # truncate field unless otherwise specified
        if truncate
            idx = sim.dis.X_idx
        else
            idx = 1:prod(sim.dis.N)
        end
    else
        idx = Int[]
    end

    # 1d plot or 2d plot
    if 1 ∈ sim.dis.N
        (sim, ψ[idx,:], theme, by, seriestype, 1)
    else
        (sim, ψ[idx,:], theme, by, seriestype, 2, 2)
    end
end


# 2d plot
@recipe function f(sim::Simulation, ψ::Array{ComplexF64}, theme::Symbol, by::Union{Function,Nothing}, seriestype::Symbol, dim1::Int, dim2::Int)

    # check whether ψ was originally provided or if just plotting sim
    if isempty(ψ)
        N=0
    else
        N = size(ψ,2)
    end

    # determine whether ψ was truncated or not and if it's just sim
    if size(ψ,1) == prod(sim.dis.N) || iszero(N)
        M = sim.dis.N
        x = sim.dis.x[1][:]
        y = sim.dis.x[2][:]
        ∂Ω = sim.bnd.∂Ω
        ε = sim.sys.ε
        F = sim.sys.F
    else
        M = sim.dis.N_tr
        x = sim.dis.x_tr[1][:]
        y = sim.dis.x_tr[2][:]
        ∂Ω = sim.bnd.∂Ω_tr
        ε = reshape(sim.sys.ε[sim.dis.X_idx],M[1],M[2])
        F = reshape(sim.sys.F[sim.dis.X_idx],M[1],M[2])
    end

    ψ_plot = Array{ComplexF64}(undef, M[1], M[2], N)
    for i in 1:N
        ψ_plot[:,:,i] = reshape(ψ[:,i], M[1], M[2])
    end

    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(theme)

    aspect_ratio :=1
    xlim := [∂Ω[1],∂Ω[2]]
    ylim := [∂Ω[3],∂Ω[4]]
    seriestype --> seriestype
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 6

    if by==nothing || iszero(N)

        # construct size
        tmp = [3,(1+N)*(∂Ω[2,2]-∂Ω[1,2])/(∂Ω[2,1]-∂Ω[1,1])]
        max_ind = findmax(tmp)[2]
        if max_ind==2
            b = 850
            a = b*tmp[1]/tmp[2]
        else
            a = 1400
            b = a*tmp[2]/tmp[1]
        end
        size --> (a,b)
        layout --> (1+N,3)

        # plot simulation
        @series begin
            title := "real"
            seriestype := :heatmap
            subplot := 1
            color := cmapsim1
            n₁ = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
            clims := (minimum(n₁), n_mult*maximum(n₁))
            x, y, permutedims(n₁)
        end
        @series begin
            title := "imag"
            seriestype := :heatmap
            subplot := 2
            color := cmapsim2
            n₂ = -F_sign*imag(sqrt.(ɛ-1im*sim.tls.D₀*F))
            clims := (-maximum(abs.(n₂)), maximum(abs.(n₂)))
            x, y, permutedims(n₂)
        end
        @series begin
            title := "F/abs²"
            seriestype := :heatmap
            subplot := 3
            color := cmapsim2
            clims := (-maximum(abs.(F)), maximum(abs.(F)))
            x, y, permutedims(F_sign*F)
        end

        for i ∈ 1:N
            ψ = ψ_plot[:,:,i]
            @series begin
                subplot := 3i+1
                clims --> (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                color := cmapc
                x, y, permutedims(real(ψ))
            end
            @series begin
                subplot := 3i+2
                clims --> (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                color := cmapc
                x, y, permutedims(imag(ψ))
            end
            @series begin
                subplot := 3i+3
                clims --> (0,maximum(abs2.(ψ)))
                color := cmapk
                x, y, permutedims(abs2.(ψ))
            end
        end
    else
        if by ∈ [abs, abs2]
            cmap = cmapk
            clims --> (0,1)
        else
            cmap = cmapc
            clims --> (-1,1)
        end
        if round(Int,N/3)==N/3
            N_col=3
            N_row = ceil(Int,N/3)
        elseif round(Int,N/2)==N/2
            N_col=2
            N_row = ceil(Int,N/2)
        elseif N==1
            N_col=1
            N_row=1
        else
            N_col = 3
            N_row = ceil(Int,N/3)
        end
        layout --> (N_row,N_col)
        tmp = 250*[N_col,(N_row)*(∂Ω[2,2]-∂Ω[1,2])/(∂Ω[2,1]-∂Ω[1,1])]
        max_ind = findmax(tmp)[2]
        if max_ind==2
            b = 850
            a = b*tmp[1]/tmp[2]
        else
            a = 1400
            b = a*tmp[2]/tmp[1]
        end
        size --> (a,b)
        for i ∈ 1:N
            ψ = ψ_plot[:,:,i]
            n₁ = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
            renorm = (maximum(abs.(ψ)) - minimum(abs.(ψ)))/(maximum(n₁)-minimum(n₁))

            n₁ = n₁ .- minimum(n₁)
            n₁ = renorm*n₁
            @series begin
                subplot := i
                seriestype := :heatmap
                seriesalpha := 1.
                if by ∈ [abs, abs2]
                    color := cmapsim1
                else
                    n₁ = n₁ .- maximum(abs.(ψ))
                    color := cmapsim1
                end
                x, y, permutedims(n₁)
            end
            @series begin
                subplot := i
                color := cmap
                seriesalpha --> .85
                if by ∈ [abs, abs2]
                    clims --> (0, maximum(by.(ψ)))
                else
                    clims --> (-maximum(abs.(ψ)),+maximum(abs.(ψ)))
                end
                x, y, permutedims(by.(ψ))
            end
        end
    end
end


# 1d plot
@recipe function f(sim::Simulation, ψ::Array{ComplexF64}, theme::Symbol, by::Union{Function,Nothing}, seriestype::Symbol, dim1::Int)

    if isempty(ψ)
        N=0
    else
        N = size(ψ,2)
    end

    if size(ψ,1) == prod(sim.dis.N) || isempty(ψ)
        if sim.dis.N[1]==1
            x = sim.dis.x[2][:]
            M = sim.dis.N[2]
            ∂Ω = sim.bnd.∂Ω[:,2]
            ε = sim.sys.ε[1,:]
            F = sim.sys.F[1,:]
            bottom_label = "y"
        else
            x = sim.dis.x[1][:]
            M = sim.dis.N[1]
            ∂Ω = sim.bnd.∂Ω[:,1]
            ε = sim.sys.ε[:,1]
            F = sim.sys.F[:,1]
            bottom_label = "x"
        end
    else
        if sim.dis.N[1]==1
            x = sim.dis.x_tr[2][:]
            M = sim.dis.N_tr[2]
            ∂Ω = sim.bnd.∂Ω_tr[:,2]
            ε = sim.sys.ε[1,sim.dis.X_idx]
            F = sim.sys.F[1,sim.dis.X_idx]
            bottom_label = "y"
        else
            x = sim.dis.x_tr[1][:]
            M = sim.dis.N_tr[1]
            ∂Ω = sim.bnd.∂Ω_tr[:,1]
            ε = sim.sys.ε[sim.dis.X_idx,1]
            F = sim.sys.F[sim.dis.X_idx,1]
            bottom_label = "x"
        end
    end

    ψ_plot = ψ
    n_mult = 1.1
    seriestype := :line
    overwrite_figure := false
    legend --> false
    lw --> 2

    if by==nothing

        size --> 250*[3,1+N]
        layout := (1+N,3)

        @series begin
            ylabel := string("Re{n(", bottom_label, ")}")
            xlabel --> bottom_label
            subplot := 1
            n = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
            ylims := (minimum(n), n_mult*maximum(n))
            x, n
        end
        @series begin
            ylabel := string("Im{n(", bottom_label, ")}")
            xlabel --> bottom_label
            subplot := 2
            n = imag(sqrt.(ɛ-1im*sim.tls.D₀*F))
            ylims := (-maximum(abs.(n)), maximum(abs.(n)))
            x, n
        end
        @series begin
            ylabel := string("F(", bottom_label, ")")
            xlabel --> bottom_label
            subplot := 3
            ylims := (-maximum(abs.(F)), maximum(abs.(F)))
            x, F
        end

        for i ∈ 1:N
            ψ = ψ_plot[:,i]
            @series begin
                ylabel := LaTeXString("Re\\{\\psi\\}")
                xlabel --> bottom_label
                subplot := 3i+1
                ylims = (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                x, real(ψ)
            end
            @series begin
                ylabel := LaTeXString("Im\\{\\psi\\}")
                xlabel --> bottom_label
                subplot := 3i+2
                ylims = (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                x, imag(ψ)
            end
            @series begin
                ylabel := LaTeXString("|\\psi|^2")
                xlabel --> bottom_label
                subplot := 3i+3
                ylims := (0, maximum(abs2.(ψ)))
                x, abs2.(ψ)
            end
        end
    else
        layout := N
        for i ∈ 1:N
            ψ = ψ_plot[:,i]
            @series begin
                if by ∈ [abs, abs2]
                    ylims := (0, by.(ψ))
                else
                    ylims := (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                end
                x, by.(ψ)
            end
        end
    end
end


################################################################################
########## ANIMATION
################################################################################
"""
    iterator = wave(sim, ψ, theme=:default; by=real, n=60, seriestype=:heatmap)

input for Plots.animate

Use cases:

`animate(wave(sim,ψ), file_name)` creates a .gif with filename

`animate(wave(sim,ψ; n=20), file_name, fps=10)` creates a 2 second movie

Note: default `fps`=20, and `n`=60, so default movie is 3 seconds long
"""
function wave(sim::Simulation, ψ, theme::Symbol=DEFAULT_COLOR_SCHEME::Symbol; truncate=true, by=real, n=60, seriestype=:heatmap)
    if truncate
        idx = sim.dis.X_idx
    else
        idx = 1:prod(sim.dis.N)
    end
    if 1 ∈ sim.dis.N
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ[idx,:], theme, by, seriestype, 1), 0:2π/n:2π*(1-1/n))
    else
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ[idx,:], theme, by, seriestype, 1, 1), 0:2π/n:2π*(1-1/n))
    end
end


################################################################################
########## AUXILLIARIES
################################################################################
"""
    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(theme)
"""
function fix_colormap(theme)

    F_sign = +1
    n_mult=1

    if theme ∈ [:dark, :juno]
        cmapc=:bkr
        cmapk=:ice
        cmapsim1=:dimgray
        cmapsim2=:bky
        n_mult=1.3
    elseif theme ∈ [:solarized]
        cmapc=:bky
        cmapk=:solar
        cmapsim1=:viridis
        cmapsim2=:bky
        n_mult=1.3
    elseif theme ∈ [:orange]
        cmapc=:bky
        cmapk=:haline
        cmapsim1=:inferno
        cmapsim2=:bky
        n_mult=1.1
    elseif theme ∈ [:lime]
        cmapc=:bky
        cmapk=:haline
        cmapsim1=:inferno
        cmapsim2=:bky
        n_mult=1.0
    elseif theme ∈ [:solarized_light]
        cmapc=:Spectral
        cmapk=:Greys
        cmapsim1=:Greys
        cmapsim2=:RdGy
        F_sign = -1
    else
        cmapc=:RdBu
        cmapk=:Greys
        cmapsim1=:Greys
        cmapsim2=:RdGy
        F_sign = -1
    end

    return cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign
end
