#TODO: check 1d plots, to waveguixe_dispersion
#TODO: make things depend on environmental variable COLOR_SCHEME

################################################################################
########## SIMULATION
################################################################################
"""
    p = plot(sim)
"""
@recipe function f(sim::Simulation)
    println("here")
    (sim, ComplexF64[])
end


################################################################################
########## SOLUTIONS
################################################################################
"""
    p = plot(sim, ψ; by=nothing, truncate=true)

plots
to turn off translucent effect, add optional argument `seriesalpha=0`
vary type of plot with `seriestype`, e.g. `seriestype=:contour`
"""
@recipe function f(sim::Simulation, ψ::Array; by=nothing, truncate=true)
println("here")
    # check whether ψ was originally provided or if just plotting sim
    if !isempty(ψ)
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
        (sim, ψ[idx,:], by, 1)
    else
        (sim, ψ[idx,:], by, 2, 2)
    end
end



# 1d plot
@recipe function f(sim::Simulation, ψ::Array, by::Union{Function,Nothing}, dim1::Int)

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



# 2d plot
@recipe function f(sim::Simulation, ψ::Array{ComplexF64}, by::Union{Function,Nothing}, dim1::Int, dim2::Int)

    # check whether ψ was originally provided or if just plotting sim
    if isempty(ψ)
        N=0
    else
        N = size(ψ,2)
    end

    # determine whether ψ was truncated or not and if it's just sim
    if size(ψ,1) == prod(sim.dis.N) || iszero(N)
        M = sim.dis.N
        x, y = sim.dis.x[1][:], sim.dis.x[2][:]
        ∂Ω = sim.bnd.∂Ω
        ε, F = sim.sys.ε, sim.sys.F
    else
        M = sim.dis.N_tr
        x, y = sim.dis.x_tr[1][:], sim.dis.x_tr[2][:]
        ∂Ω = sim.bnd.∂Ω_tr
        ε = reshape(sim.sys.ε[sim.dis.X_idx],M[1],M[2])
        F = reshape(sim.sys.F[sim.dis.X_idx],M[1],M[2])
    end

    ψ_plot = Array{ComplexF64}(undef, M[1], M[2], N)
    for i in 1:N
        ψ_plot[:,:,i] = reshape(ψ[:,i], M[1], M[2])
    end

    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(:default)

    aspect_ratio :=1
    xlim := [∂Ω[1],∂Ω[2]]
    ylim := [∂Ω[3],∂Ω[4]]
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 6
    seriestype --> :heatmap
    legend --> false

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
        for j ∈ 1:3
            for i ∈ CartesianIndices(sim.bnd.bc)
                if sim.bnd.bc[i] == :d
                    linestyle := :solid
                elseif sim.bnd.bc[i] == :n
                    linestype := :dash
                end
                if i[2]==1
                    data = ([sim.bnd.∂Ω[i], sim.bnd.∂Ω[i]], [sim.bnd.∂Ω[1,2], sim.bnd.∂Ω[2,2]])
                else
                    data = ([sim.bnd.∂Ω[1,1], sim.bnd.∂Ω[2,1]], [sim.bnd.∂Ω[i], sim.bnd.∂Ω[i]])
                end
                if sim.bnd.bl[i] !== :none
                    if sim.bnd.bl[i] ∈ [:pml_out, :abs_out]
                        fill_color = :red
                    elseif sim.bnd.bl[i] ∈ [:pml_in, :abs_in]
                        fill_color = :blue
                    end
                    if i[2] == 1
                        datax = [ sim.bnd.∂Ω[i[1],1], sim.bnd.∂Ω_tr[i[1],1], sim.bnd.∂Ω_tr[i[1],1], sim.bnd.∂Ω[i[1],1], sim.bnd.∂Ω[i[1],1] ]
                        datay = [ sim.bnd.∂Ω[1,2], sim.bnd.∂Ω[1,2],    sim.bnd.∂Ω[2,2],    sim.bnd.∂Ω[2,2], sim.bnd.∂Ω[1,2] ]
                    else
                        datax = [ sim.bnd.∂Ω[1,1], sim.bnd.∂Ω[1,1], sim.bnd.∂Ω[2,1], sim.bnd.∂Ω[2,1], sim.bnd.∂Ω[1,1] ]
                        datay = [ sim.bnd.∂Ω[i[1],2], sim.bnd.∂Ω_tr[i[1],2], sim.bnd.∂Ω_tr[i[1],2], sim.bnd.∂Ω[i[1],2], sim.bnd.∂Ω[i[1],2] ]
                    end
                    @series begin
                        seriestype := :path
                        subplot := j
                        lw := 0
                        fill := (0,.15,fill_color)
                        (datax, datay)
                    end
                end
                if sim.bnd.bc[i] ∈ [:d, :n]
                    @series begin
                        seriestype := :line
                        color := :black
                        subplot := j
                        lw := 2
                        data
                    end
                end
            end
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





################################################################################
########## ANIMATION
################################################################################
"""
    iterator = wave(sim, ψ; by=real, n=60, seriestype=:heatmap)

input for Plots.animate

Use cases:

`animate(wave(sim,ψ), file_name)` creates a .gif with filename

`animate(wave(sim,ψ; n=20), file_name, fps=10)` creates a 2 second movie

Note: default `fps`=20, and `n`=60, so default movie is 3 seconds long
"""
function wave(sim::Simulation, ψ; truncate=true, by=real, n=60)
    if truncate
        idx = sim.dis.X_idx
    else
        idx = 1:prod(sim.dis.N)
    end
    if 1 ∈ sim.dis.N
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ[idx,:], by, 1), 0:2π/n:2π*(1-1/n))
    else
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ[idx,:], by, 1, 1), 0:2π/n:2π*(1-1/n))
    end
end
