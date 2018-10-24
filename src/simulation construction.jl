#    TODO: sub_pixel_smoothing in 2-dim case, take care of smoothing along edges
#    TODO: fix Bravais in 1-dimension
    #TODO: build_channels! higher order band gaps
################################################################################
###############   SUB PIXEL SMOOTHING    #######################################
################################################################################
"""
    ε, F, regions = sub_pixel_smoothing(bnd, dis, sys)
"""
function sub_pixel_smoothing(bnd::Boundary, dis::Discretization, sys::System; disp_opt=false)

    x = dis.x[1]
    y = dis.x[2]

    domains = ScalarFDFD.which_domain.(x, y, Ref(bnd), Ref(sys))
    xy = ScalarFDFD.bravais_coordinates_unit_cell.(x, y, domains, Ref(sys))
    xb = Array{Float64}(undef,size(xy))
    yb = Array{Float64}(undef,size(xy))
    for i ∈ eachindex(xb)
        xb[i] = xy[i][1]
        yb[i] = xy[i][2]
    end
    regions = ScalarFDFD.which_region.(xb, yb, domains, Ref(sys))
    r = regions

    ε_by_region = sys.ε_by_region
    F_by_region = sys.F_by_region

    ɛ = ɛ_by_region[regions]
    F = F_by_region[regions]

    sub_pixel_num = dis.sub_pixel_num
    if length(y) == 1
        sub_x = Array{Float64}(undef, sub_pixel_num,1)
        sub_y = y
        xy = Array{Tuple{Float64,Float64}}(undef, sub_pixel_num, 1)
        xb = Array{Float64}(undef, sub_pixel_num, 1)
        yb = Array{Float64}(undef, sub_pixel_num, 1)
        sub_regions = Array{Int}(undef, sub_pixel_num, 1)
        sub_domains = Array{Int}(undef, sub_pixel_num, 1)
        # pg = Progress((length(x)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
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
                F[i] = mean(F_by_region[sub_regions])
            end
            # next!(pg)
        end
    elseif length(x) == 1
        sub_x = x
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xy = Array{Tuple{Float64,Float64}}(undef, 1, sub_pixel_num)
        xb = Array{Float64}(undef, 1, sub_pixel_num)
        yb = Array{Float64}(undef, 1, sub_pixel_num)
        sub_regions = Array{Int}(undef, 1, sub_pixel_num)
        sub_domains = Array{Int}(undef, 1, sub_pixel_num)
        # pg = Progress((length(y)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
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
            # next!(pg)
        end
    else
        # rng = MersenneTwister(0)
        # Rx = rand(rng,sub_pixel_num)
        # Ry = rand(rng,sub_pixel_num)
        sub_x = Array{Float64}(undef, sub_pixel_num, 1)
        sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        xy = Array{Tuple{Float64,Float64}}(undef, sub_pixel_num, sub_pixel_num)
        xb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        yb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        sub_regions = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        sub_domains = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
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
                # Rx = rand(rng,sub_pixel_num)
                # Ry = rand(rng,sub_pixel_num)
                sub_x[:] = LinRange(x_min, x_max, sub_pixel_num) #x_min .+ (x_max-x_min)*Rx
                sub_y[:] = LinRange(y_min, y_max, sub_pixel_num) #y_min .+ (y_max-y_min)*Ry
                sub_domains[:] = ScalarFDFD.which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
                xy[:] = ScalarFDFD.bravais_coordinates_unit_cell.(sub_x, sub_y, sub_domains, Ref(sys))
                for k ∈ eachindex(xb)
                    xb[k] = xy[k][1]
                    yb[k] = xy[k][2]
                end
                sub_regions[:] = ScalarFDFD.which_region.(xb, yb, sub_domains, Ref(sys))
                ɛ[i,j] = mean(ɛ_by_region[sub_regions])
                F[i,j] = mean(F_by_region[sub_regions])
            end
            if disp_opt
                next!(pg)
            end
        end
    end

    return ɛ, F, regions
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
    while i≤length(sys.domains[domain].is_in_region) && !sys.domains[domain].is_in_region[i](x, y, i, sys.domains[domain])
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
function whole_domain(x,y,idx,bnd::Boundary,sys::System)
    return true
end
