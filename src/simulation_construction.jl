# TODO: sub_pixel_smoothing in 2-dim case, take care of smoothing along edges
# TODO: fix Bravais in 1-dimension

################################################################################
###############   SUB PIXEL SMOOTHING    #######################################
################################################################################
"""
    ε, regions = sub_pixel_smoothing(bnd, dis, sys; disp_opt=false)
"""
function sub_pixel_smoothing(bnd::Boundary, dis::Discretization, sys::System; disp_opt=false)

    x = dis.X
    y = dis.Y

    domains = which_domain.(x, y, Ref(bnd), Ref(sys))
    xb, yb = bravais_coordinates_unit_cell(x, y, domains, sys)
    r = which_region.(xb, yb, domains, Ref(sys))

    ε = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
    for i ∈ CartesianIndices(r)
        ε[i] = sys.ε_by_region[r[i]](x[i[1]],y[i[2]],sys.params_by_region[r[i]])
    end

    sub_pixel_num = dis.sub_pixel_num
    if dis.N[2] == 1
        sub_x = Array{Float64}(undef, sub_pixel_num,1)
        sub_y = y
        xb = Array{Float64}(undef, sub_pixel_num, 1)
        yb = Array{Float64}(undef, sub_pixel_num, 1)
        sub_regions = Array{Int}(undef, sub_pixel_num, 1)
        sub_domains = Array{Int}(undef, sub_pixel_num, 1)
        for i ∈ 2:(size(x,1)-1)
            nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
            if nearestNeighborFlag
                x_min = (x[i]+x[i-1])/2
                x_max = (x[i]+x[i+1])/2
                sub_x[:] = LinRange(x_min, x_max, sub_pixel_num)
                sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, sub_domains, Ref(sys))
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                ɛ[i] = mean(ɛ_by_region[sub_regions])
            end
        end
    elseif dis.N[1] == 1
        sub_x = x
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xb = Array{Float64}(undef, 1, sub_pixel_num)
        yb = Array{Float64}(undef, 1, sub_pixel_num)
        sub_regions = Array{Int}(undef, 1, sub_pixel_num)
        sub_domains = Array{Int}(undef, 1, sub_pixel_num)
        for i ∈ 2:(size(y,2)-1)
            nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
            if nearestNeighborFlag
                y_min = (y[i]+y[i-1])/2
                y_max = (y[i]+y[i+1])/2
                sub_y[:] = LinRange(y_min, y_max, sub_pixel_num)
                sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, sub_domains, Ref(sys))
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                ɛ[i] = mean(ɛ_by_region[sub_regions])
                F[i] = mean(F_by_region[sub_regions])
            end
        end
    else
        sub_x = Array{Float64}(undef, sub_pixel_num, 1)
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        yb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        sub_regions = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        sub_domains = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        ε_temp = Array{ComplexF64}(undef,sub_pixel_num, sub_pixel_num)
        if disp_opt
            pg = Progress((size(x,1)-2)*(size(y,2)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
        end
        for i ∈ 2:(size(x,1)-1), j ∈ 2:(size(y,2)-1)
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
                bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, sub_domains, Ref(sys))
                sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
                for i ∈ CartesianIndices(sub_regions)
                    ε_temp[i] = sys.ε_by_region[r[i]](x[i[1]],y[i[2]],sys.params_by_region[r[i]])
                end
                ε[i,j] = mean(ε_temp)
            end
            if disp_opt
                next!(pg)
            end
        end
    end

    return ɛ, r
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
    @assert domain≤length(sys.domains) "no domain found for ($x,$y). has a background domain been defined?"
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
        if isDirichlet(bc[i])
            bc[i] = :d
        elseif isNeumann(bc[i])
            bc[i] = :n
        elseif isOpen(bc[i])
            bc[i] = :o
        elseif isPeriodic(bc[i])
            bc[i] = :p
        else
            throw("invalid boundary condition $(bc[i])")
        end
    end
    return nothing
end


"""
    fix_bl!(bl)
"""
function fix_bl!(bl)
    for i ∈ eachindex(bl)
        if isPMLout(bl[i])
            bl[i] = :pml_out
        elseif isPMLin(bl[i])
            bl[i] = :pml_in
        elseif isABSout(bl[i])
            bl[i] = :abs_out
        elseif isABSin(bl[i])
            bl[i] = :abs_in
        elseif isNone(bl[i])
            bl[i] = :none
        else
            throw("invalid boundary layer $(bl[i])")
        end
    end
    return nothing
end


################################################################################
###############  STANDARD DOMAINS  #######################################
################################################################################
"""
    left_domain(x, y, domain_index, bnd, sys)
"""
function left_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :left && x<bnd.∂Ω_tr[1,1]
end


"""
    right_domain(x, y, domain_index, bnd, sys)
"""
function right_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :right && x>bnd.∂Ω_tr[2,1]
end


"""
    bottom_domain(x, y, domain_index, bnd, sys)
"""
function bottom_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :bottom && y<bnd.∂Ω_tr[1,2]
end


"""
    top_domain(x, y, domain_index, bnd, sys)
"""
function top_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :top && y>bnd.∂Ω_tr[2,2]
end


"""
    whole_domain(x, y, domain_index, bnd, sys)
"""
function whole_domain(x,y,idx::Int,bnd::Boundary,sys::System)::Bool
    return true
end


"""
    piecewise_constant_ε(x,y,params)
"""
function piecewise_constant_ε(x,y,params)
    return complex(params[:n₁],params[:n₂])^2
end



################################################################################################
### BOUNDARY LAYERS
################################################################################################
"""
    Σd, Σe = σ(x, y, bnd, side::Symbol)
    Σd, Σe = σ(dis, bnd)
    Σd, Σe = σ(sim)

conductivity for absorbing layer (PML or not). `1/(1+Σd/k)` is sandwiched between derivatives in laplacian, `(1+Σe/k) multiplies ε`
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
         u = abs(pos[j]-bnd.∂Ω_tr[i,j])/(bnd.bl_depth[i,j] + eps())
         β = 2*sqrt(bnd.bl_depth[i,j]) + eps()
         Σ = α[i,j]*u*expm1(β*u)/expm1(β)
     else
         Σ = complex(0.)
    end
    return Σ
end
function σ(dis::Discretization, bnd::Boundary)
    Σ1 = zeros(ComplexF64,dis.N[1]+1,dis.N[2]+1)
    Σ2 = zeros(ComplexF64,dis.N[1]  ,dis.N[2]  )
    x1 = vcat(dis.x[1] .- dis.dx[1]/2, dis.x[1][end] + dis.dx[1]/2)
    y1 = hcat(dis.x[2] .- dis.dx[2]/2, dis.x[2][end] + dis.dx[2]/2)
    x2, y2 = dis.x[1], dis.x[2]
    if isPolar(dis.coordinate_system)
        X1, Y1 = x1.*cos.(y1), x1.*sin.(y1)
        X2, Y2 = x2.*cos.(y2), x2.*sin.(y2)
    else
        X1, Y1 = broadcast((x,y)->x, x1,y1), broadcast((x,y)->y, x1,y1)
        X2, Y2 = broadcast((x,y)->x, x2,y2), broadcast((x,y)->y, x2,y2)
    end
    sides = [:left, :right, :bottom, :top]
    for s ∈ sides
        Σ1 += σ.(X1, Y1, Ref(bnd), s)
        Σ2 += σ.(X2, Y2, Ref(bnd), s)
    end
    return Σ1, Σ2
end
function σ(sim::Simulation)
    return σ(sim.dis, sim.bnd)
end
