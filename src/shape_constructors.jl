#TODO: extend signature of circle to all other shapes (i.e. add to sim, sys, dom, etc)
# """
#     1. simulation = add_circle(sim; R, x0, y0, n₁=1, n₂=0, F=0)
#
#     2. system = add_circle(sys; R, x0, y0, n₁=1, n₂=0, F=0)
#
#     3. domain = add_circle(dom; R, x0, y0, n₁=1, n₂=0, F=0)
#
#     4. regions = add_circle(old_regions=Regions(); R, x0, y0, n₁=1, n₂=0, F=0)
#
# 1. add a circular domain to a Simulation
# 2. add a circular domain to a System
# 3. add a circular region to a Domain
# 4. add a circular region to a Region object,
#
# all with index and pump specified by `n₁`, etc.
# """

################################################################################
### CIRCLE
################################################################################
"""
    new_domain = add_circle_region(domain; R, x0, y0, n₁=1, n₂=0, F=0))

    new_regions = add_circle_region(old_regions=Regions(); R, x0, y0, n₁=1, n₂=0, F=0))
"""
function add_circle_region(domain::Domain; R, x0, y0, n₁=1, n₂=0, F=0)
    return add_region(circle_region, [x0, y0, R], domain, n₁, n₂, F)
end
function add_circle_region(old_regions::Regions=Regions(); R, x0, y0, n₁=1, n₂=0, F=0)
    return add_region(circle_region, [x0, y0, R], n₁, n₂, F, old_regions)
end

"""
    new_domain = circle_domain(; R, x0, y0, n₁=1, n₂=0, F=0)
"""
function circle_domain(; R, x0, y0, n₁=1, n₂=0, F=0)
    return build_domain(Regions(), is_in_domain=circle_region, domain_params=[x0, y0, R], n₁=n₁, n₂=n₂, F=F, domain_type=:circle)
end

"""
    modified_sim = add_circle_domain(sim, R, x0, y0, n₁=1, n₂=0, F=0)

    modified_sys = add_circle_domain(sys, R, x0, y0, n₁=1, n₂=0, F=0)
"""
function add_circle_domain(sim::Simulation; R, x0, y0, n₁=1, n₂=0, F=0)
    sys = circle_domain(sim.sys; R=R, x0=x0, y0=y0, n₁=n₁, n₂=n₂, F=F)
    return Simulation(sim; :sys => sys)
end
function add_circle_domain(sys::System; R, x0, y0, n₁=1, n₂=0, F=0)
    domain = circle_domain(Regions(), is_in_domain=circle_region, domain_params=[x0, y0, R], n₁=n₁, n₂=n₂, F=F, domain_type=:circle)
    return System(vcat(domain,sys.domains))
end



################################################################################
### RECTANGLE
################################################################################
"""
    new_regions = add_rectangle_region(old_regions=Regions(); a, b, x0, y0, θ, n₁=1, n₂=0, F=0)

    modified_domain = add_rectangle_region(domain; a, b, x0, y0, θ, n₁=1, n₂=0, F=0)
"""
function add_rectangle_region(domain::Domain; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
    return add_region(rectangle_region, [x0, y0, cos(α), sin(α), a, b], domain, n₁, n₂, F, old_regions)
end
function add_rectangle_region(old_regions::Regions=Regions(); a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
    return add_region(rectangle_region, [x0, y0, cos(α), sin(α), a, b], n₁, n₂, F, old_regions)
end
"""
    new_domain = rectangle_domain(; a, b, x0, y0, α=0, , n₁=1, n₂=0, F=0)
"""
function rectangle_domain(; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
    return build_domain(Regions(), is_in_domain=rectangle_region, domain_params=[x0, y0, cos(α), sin(α), a, b], n₁=n₁, n₂=n₂, F=F, domain_type=:rectangle)
end
"""
    modified_sim = add_rectangle_domain(sim; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)

    modified_sys = add_rectangle_domain(sys; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
"""
function add_rectangle_domain(sim::Simulation; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
    sys = rectangle_domain(sim.sys; a=a, b=b, x0=x0, y0=y0, α=α, n₁=n₁, n₂=n₂, F=F)
    return Simulation(sim; :sys => sys)
end
function add_rectangle_domain(sys::System; a, b, x0, y0, α=0, n₁=1, n₂=0, F=0)
    domain = rectangle_domain(; a=a, b=b, x0=x0, y0=y0, α=α, n₁=n₁, n₂=n₂, F=F)
    return System(vcat(domain,sys.domains))
end




"""
    new_regions = add_ellipse(old_regions=Regions(); a, b, x0, y0, θ, n₁=1, n₂=0, F=0)

    new_domain = add_ellipse(domain; a, b, x0, y0, θ, n₁=1, n₂=0, F=0)
"""
function add_ellipse(old_regions::Regions=Regions(); a, b, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(ellipse_region, [x0, y0, cos(α), sin(α), a, b], n₁, n₂, F, old_regions)
end
function add_ellipse(domain::Domain; a, b, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(ellipse_region, [x0, y0, cos(α), sin(α), a, b], domain, n₁, n₂, F)
end


"""
    new_regions = add_square(old_regions=Regions(); a, x0, y0, θ, n₁=1, n₂=0, F=0)

    new_domain = add_square(domain; a, x0, y0, θ, n₁=1, n₂=0, F=0)
"""
function add_square(old_regions::Regions=Regions(); a, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(square_region, [x0, y0, cos(α), sin(α), a], n₁, n₂, F, old_regions)
end
function add_square(domain::Domain; a, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(square_region, [x0, y0, cos(α), sin(α), a], domain, n₁, n₂, F)
end



"""
    new_regions = add_parallelogram(old_regions=Regions(); a, b, α, x0, y0, θ, n₁=1, n₂=0, F=0)

    new_domain = add_parallelogram(domain; a, b, α, x0, y0, θ, n₁=1, n₂=0, F=0)
"""
function add_parallelogram(old_regions::Regions=Regions(); a, b, θ, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(parallelogram_region, [x0, y0, cos(α), sin(α), a, b, cos(θ), tan(θ)], n₁, n₂, F, old_regions)
end
function add_parallelogram(domain::Domain; a, b, θ, x0, y0, α, n₁=1, n₂=0, F=0)
    return add_region(parallelogram_region, [x0, y0, cos(α), sin(α), a, b, cos(θ), tan(θ)], domain, n₁, n₂, F)
end
