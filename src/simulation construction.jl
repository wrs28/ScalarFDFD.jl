#    TODO: sub_pixel_smoothing in 2-dim case, take care of smoothing along edges
#    TODO: fix Bravais in 1-dimension
#TODO: build_channels! higher order band gaps
################################################################################
###############   SUB PIXEL SMOOTHING    #######################################
################################################################################

"""
    ε, regions = sub_pixel_smoothing(bnd, dis, sys; disp_opt=false)
"""
function sub_pixel_smoothing(bnd::Boundary, dis::Discretization, sys::System; disp_opt=false)

    x = dis.x[1]
    y = dis.x[2]

    domains = which_domain.(x, y, Ref(bnd), Ref(sys))
    xy = bravais_coordinates_unit_cell.(x, y, domains, Ref(sys))
    xb = Array{Float64}(undef,size(xy))
    yb = Array{Float64}(undef,size(xy))
    for i ∈ eachindex(xb)
        xb[i] = xy[i][1]
        yb[i] = xy[i][2]
    end
    regions = which_region.(xb, yb, domains, Ref(sys))
    r = regions

    ε = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
    for i ∈ CartesianIndices(regions)
        ε[i] = sys.ε_by_region[regions[i]](x[i[1]],y[i[2]],sys.params_by_region[regions[i]])
    end

    sub_pixel_num = dis.sub_pixel_num
    if length(y) == 1
        sub_x = Array{Float64}(undef, sub_pixel_num,1)
        sub_y = y
        xy = Array{Tuple{Float64,Float64}}(undef, sub_pixel_num, 1)
        xb = Array{Float64}(undef, sub_pixel_num, 1)
        yb = Array{Float64}(undef, sub_pixel_num, 1)
        sub_regions = Array{Int}(undef, sub_pixel_num, 1)
        sub_domains = Array{Int}(undef, sub_pixel_num, 1)
        for i in 2:(length(x)-1)
            nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
            if nearestNeighborFlag
                x_min = (x[i]+x[i-1])/2
                x_max = (x[i]+x[i+1])/2
                sub_x[:] = LinRange(x_min, x_max, sub_pixel_num)
                sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                xy[:] = bravais_coordinates_unit_cell.(sub_x, sub_y, sub_domains, Ref(sys))
                for k ∈ eachindex(xb)
                    xb[k] = xy[k][1]
                    yb[k] = xy[k][2]
                end
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                ɛ[i] = mean(ɛ_by_region[sub_regions])
            end
        end
    elseif length(x) == 1
        sub_x = x
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xy = Array{Tuple{Float64,Float64}}(undef, 1, sub_pixel_num)
        xb = Array{Float64}(undef, 1, sub_pixel_num)
        yb = Array{Float64}(undef, 1, sub_pixel_num)
        sub_regions = Array{Int}(undef, 1, sub_pixel_num)
        sub_domains = Array{Int}(undef, 1, sub_pixel_num)
        for i in 2:(length(y)-1)
            nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
            if nearestNeighborFlag
                y_min = (y[i]+y[i-1])/2
                y_max = (y[i]+y[i+1])/2
                sub_y[:] = LinRange(y_min, y_max, sub_pixel_num)
                sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                xy[:] = bravais_coordinates_unit_cell.(sub_x, sub_y, sub_domains, Ref(sys))
                for k ∈ eachindex(xb)
                    xb[k] = xy[k][1]
                    yb[k] = xy[k][2]
                end
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                ɛ[i] = mean(ɛ_by_region[sub_regions])
                F[i] = mean(F_by_region[sub_regions])
            end
        end
    else
        sub_x = Array{Float64}(undef, sub_pixel_num, 1)
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xy = Array{Tuple{Float64,Float64}}(undef, sub_pixel_num, sub_pixel_num)
        xb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        yb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        sub_regions = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        sub_domains = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        ε_temp = Array{ComplexF64}(undef,sub_pixel_num, sub_pixel_num)
        if disp_opt
            pg = Progress((length(x)-2)*(length(y)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
        end
        for i in 2:(length(x)-1), j in 2:(length(y)-1)
            nearestNeighborFlag = r[i,j]!==r[i,j+1] || r[i,j]!==r[i,j-1] || r[i,j]!==r[i+1,j] || r[i,j]!==r[i-1,j]
            nextNearestNeighborFlag = r[i,j]!==r[i+1,j+1] || r[i,j]!==r[i-1,j-1] || r[i,j]!==r[i+1,j-1] || r[i,j]!==r[i-1,j+1]
            if nearestNeighborFlag || nextNearestNeighborFlag
                x_min = (x[i]+x[i-1])/2
                y_min = (y[j]+y[j-1])/2
                x_max = (x[i]+x[i+1])/2
                y_max = (y[j]+y[j+1])/2
                sub_x[:] = LinRange(x_min, x_max, sub_pixel_num)
                sub_y[:] = LinRange(y_min, y_max, sub_pixel_num)
                sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                xy[:] = bravais_coordinates_unit_cell.(sub_x, sub_y, sub_domains, Ref(sys))
                for k ∈ eachindex(xb)
                    xb[k] = xy[k][1]
                    yb[k] = xy[k][2]
                end
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                for i ∈ CartesianIndices(sub_regions)
                    ε_temp[i] = sys.ε_by_region[regions[i]](x[i[1]],y[i[2]],sys.params_by_region[regions[i]])
                end
                ε[i,j] = mean(ε_temp)
            end
            if disp_opt
                next!(pg)
            end
        end
    end

    return ɛ, regions
end


"""
    which_domain(x, y, bnd, sys)
"""
function which_domain(x, y, bnd::Boundary, sys::System)
    domain = 1
    while (domain ≤ length(sys.domains) &&
        !(
        (left_domain(x,y,domain,bnd,sys) && sys.domains[domain].is_in_domain(x, y, domain, bnd, sys)::Bool) ||
        (right_domain(x,y,domain,bnd,sys) && sys.domains[domain].is_in_domain(x, y, domain, bnd, sys)::Bool) ||
        (bottom_domain(x,y,domain,bnd,sys) && sys.domains[domain].is_in_domain(x, y, domain, bnd, sys)::Bool) ||
        (top_domain(x,y,domain,bnd,sys) && sys.domains[domain].is_in_domain(x, y, domain, bnd, sys)::Bool) ||
        (sys.domains[domain].which_asymptote==:none && sys.domains[domain].is_in_domain(x, y, domain, bnd, sys)::Bool)
        ))
        domain += 1
    end
    if domain == length(sys.domains)+1
        throw(ErrorException("no domain found for ($x,$y). has a background domain been defined?"))
    end
    return domain
end



"""
    region = which_region(x, y, domain_index, sys)
"""
function which_region(x, y, domain, sys::System)
    i = 1
    while i≤length(sys.domains[domain].is_in_subdomain) && !sys.domains[domain].is_in_subdomain[i](x, y, i, sys.domains[domain])
        i += 1
    end
    region = sys.num_prev_regions[domain] + i
    return region
end


################################################################################
###############  SYSTEM STANDARDIZATION  #######################################
################################################################################
"""
    fix_bc!(bc)
"""
function fix_bc!(bc)
    for i ∈ eachindex(bc)
        if bc[i] ∈ [:d, :D, :dirichlet, :Dirichlet, :hard, :h]
            bc[i] = :d
        elseif bc[i] ∈ [:n, :N, :neumann, :Neumann, :soft, :s]
            bc[i] = :n
        elseif bc[i] ∈ [:o, :open, :Open]
            bc[i] = :o
        elseif bc[i] ∈ [:p, :periodic, :Periodic, :bloch, :Bloch]
            bc[i] = :p
        else
            throw(ErrorException("invalid boundary condition $(bc[i])"))
        end
    end
    return nothing
end


"""
    fix_bl!(bl)
"""
function fix_bl!(bl)
    for i ∈ eachindex(bl)
        if bl[i] ∈ [:pml_out, :PML_OUT, :PML_out, :pml, :PML]
            bl[i] = :pml_out
        elseif bl[i] ∈ [:pml_in, :PML_IN, :PML_in, :pml_conj, :PML_conj]
            bl[i] = :pml_in
        elseif bl[i] ∈ [:abs_out, :ABS_OUT, :ABS_out, :abs, :ABS, :amp, :AMP]
            bl[i] = :abs_out
        elseif bl[i] ∈ [:abs_in, :ABS_IN, :ABS_in, :abs_conj, :ABS_conj]
            bl[i] = :abs_in
        elseif bl[i] ∈ [:none, :nothing, :empty]
            bl[i] = :none
        else
            throw(ErrorException("invalid boundary layer $(bl[i])"))
        end
    end
    return nothing
end



"""
    SA = absorbing_boundary_layers(sim, k)
"""
function absorbing_boundary_layers(sim::Simulation, k)

    Σ = sim.sys.Σ
    N = sim.dis.N
    ε = sim.sys.ε

    SA = [sparse(1:N[j], 1:N[j], Vector{ComplexF64}(undef,N[j]), N[j], N[j]) for i ∈ 1:2, j ∈ 1:2]
    for r ∈ CartesianIndices(SA)
        j = r[2]
        if sim.bnd.bl[r] ∈ [:pml_out, :pml_in]
            SA[r] = spzeros(ComplexF64,N[j],N[j])   #sparse(complex(1.,0)I, N[j], N[j])
        else
            SA[r] = sparse(1:N[j], 1:N[j], 0 .+ 1im*Σ[r]/real(k), N[j], N[j])
        end
    end

    return SA
end




################################################################################
###############  STANDARD DOMAINS  #######################################
################################################################################
"""
    left_domain(x, y, domain_index, bnd, sys)
"""
function left_domain(x,y,idx,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :left && x<bnd.∂Ω_tr[1,1]
end


"""
    right_domain(x, y, domain_index, bnd, sys)
"""
function right_domain(x,y,idx,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :right && x>bnd.∂Ω_tr[2,1]
end


"""
    bottom_domain(x, y, domain_index, bnd, sys)
"""
function bottom_domain(x,y,idx,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :bottom && y<bnd.∂Ω_tr[1,2]
end


"""
    top_domain(x, y, domain_index, bnd, sys)
"""
function top_domain(x,y,idx,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :top && y>bnd.∂Ω_tr[2,2]
end


"""
    whole_domain(x, y, domain_index, bnd, sys)
"""
function whole_domain(x,y,idx,bnd::Boundary,sys::System)::Bool
    return true
end



"""
    piecewise_constant_ε(x,y,params)
"""
function piecewise_constant_ε(x,y,params)
    return complex(params[:n₁],params[:n₂])^2
end


"""
    Σ::ComplexF64 = σ(x, y, bnd, side::Symbol)
    Σ::Array{ComplexF64,2} = σ(dis, bnd)
    Σ::Array{ComplexF64,2} = σ(sim)

conductivity for absorbing layer (PML or not)
"""
function σ(x::Real,y::Real,bnd::Boundary,side::Symbol)

    α = bnd.α
    if side == :left
        i = 1; j = 1
    elseif side == :right x > bnd.∂Ω_tr[2,1]
        i = 2; j = 1
    elseif side == :bottom
        i = 1; j = 2
    elseif side == :top
        i = 2; j = 2
    else
        throw(ArgumentError("invalid side $side, must be one of :left, :right, :bottom, :top"))
    end
    pos = (x,y)
    if sign(bnd.∂Ω_tr[i,j] - pos[j])*(-1)^i ≤ 0
         u = abs(pos[j]-bnd.∂Ω_tr[i,j])/bnd.bl_depth[i,j]
         β = 2*sqrt(bnd.bl_depth[i,j])
         Σ = α[i,j]*u*(exp(β*u)-1)/(exp(β)-1)
     else
         Σ = complex(0.)
    end
    return Σ
end
function σ(dis::Discretization, bnd::Boundary)
    Σ = zeros(ComplexF64,dis.N[1],dis.N[2])
    Σ += σ.(dis.x[1], dis.x[2], Ref(bnd), :left)
    Σ += σ.(dis.x[1], dis.x[2], Ref(bnd), :right)
    Σ += σ.(dis.x[1], dis.x[2], Ref(bnd), :bottom)
    Σ += σ.(dis.x[1], dis.x[2], Ref(bnd), :top)
    return Σ
end
function σ(sim::Simulation)
    return σ(sim.dis, sim.bnd)
end
