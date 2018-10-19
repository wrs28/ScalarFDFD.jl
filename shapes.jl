"""
    circle_region(x, y, x0, y0, R)
"""
function circle_region(x,y,x0,y0,R)
    return hypot((x-x0), (y-y0)) ≤ R
end
function circle_region(x, y, idx, dom::Domain)
    x0 = dom.region_params[idx][1]
    y0 = dom.region_params[idx][2]
    R = dom.region_params[idx][3]
    return circle_region(x,y,x0,y0,R)
end
function circle_region(x, y, idx, bnd::Boundary, sys::System)
    x0 = sys.domains[idx].domain_params[1]
    y0 = sys.domains[idx].domain_params[2]
    R = sys.domains[idx].domain_params[3]
    return circle_region(x,y,x0,y0,R)
end


"""
    ellipse_region(x, y, [a, b, x0, y0, θ])
"""
function ellipse_region(x,y,x0,y0,cosθ,sinθ,a,b)
    xrot, yrot = rotate(x-x0,y-y0,cosθ,sinθ)
    return hypot(xrot/a, yrot/b) ≤ 1
end
function ellipse_region(x, y, idx, dom::Domain)
    x0 = dom.region_params[idx][1]
    y0 = dom.region_params[idx][2]
    cosθ = dom.region_params[idx][3]
    sinθ = dom.region_params[idx][4]
    a = dom.region_params[idx][5]
    b = dom.region_params[idx][6]
    return ellipse_region(x,y,x0,y0,θ,a,b)
end
function ellipse_region(x, y, idx, bnd::Boundary, sys::System)
    x0 = sys.domains[idx].domain_params[1]
    y0 = sys.domains[idx].domain_params[2]
    cosθ = sys.domains[idx].domain_params[3]
    sinθ = sys.domains[idx].domain_params[4]
    a = sys.domains[idx].domain_params[5]
    b = sys.domains[idx].domain_params[6]
    return ellipse_region(x,y,x0,y0,θ,a,b)
end


"""
    square_region(x, y, [a, x0, y0, θ])
"""
function square_region(x,y,x0,y0,cosθ,sinθ,a)
    xrot, yrot = rotate(x-x0,y-y0,cosθ,sinθ)
    return 0 ≤ xrot ≤ a  &&  0 ≤ yrot ≤ a
end
function square_region(x, y, idx, dom::Domain)
    x0 = dom.region_params[idx][1]
    y0 = dom.region_params[idx][2]
    cosθ = dom.region_params[idx][3]
    sinθ = dom.region_params[idx][4]
    a = dom.region_params[idx][5]
    return square_region(x,y,x0,y0,cosθ,sinθ,a)
end
function square_region(x, y, idx, bnd::Boundary, sys::System)
    x0 = sys.domains[idx].domain_params[1]
    y0 = sys.domains[idx].domain_params[2]
    cosθ = sys.domains[idx].domain_params[3]
    sinθ = sys.domains[idx].domain_params[4]
    a = sys.domains[idx].domain_params[5]
    return square_region(x,y,x0,y0,cosθ,sinθ,a)
end


"""
    rectangle_region(x, y, [a, b, x0, y0, θ])
"""
function rectangle_region(x,y,x0,y0,cosθ,sinθ,a,b)
    xrot, yrot = rotate(x-x0,y-y0,cosθ,sinθ)
    return 0 ≤ xrot ≤ a  &&  0 ≤ yrot ≤ b
end
function rectangle_region(x, y, idx, dom::Domain)
    x0 = dom.region_params[idx][1]
    y0 = dom.region_params[idx][2]
    cosθ = dom.region_params[idx][3]
    sinθ = dom.region_params[idx][4]
    a = dom.region_params[idx][5]
    b = dom.region_params[idx][6]
    return rectangle_region(x,y,x0,y0,cosθ,sinθ,a,b)
end
function rectangle_region(x, y, idx, bnd::Boundary, sys::System)
    x0 = sys.domains[idx].domain_params[1]
    y0 = sys.domains[idx].domain_params[2]
    cosθ = sys.domains[idx].domain_params[3]
    sinθ = sys.domains[idx].domain_params[4]
    a = sys.domains[idx].domain_params[5]
    b = sys.domains[idx].domain_params[6]
    return rectangle_region(x,y,x0,y0,cosθ,sinθ,a,b)
end


"""
    parallelogram_region(x, y, [a, b, α, x0, y0, θ])
"""
function parallelogram_region(x,y,x0,y0,cosθ,sinθ,a,b,cosα,tanα)
    xrot, yrot = rotate(x-x0, y-y0, cosθ, sinθ)
    return tanα*(xrot-a) ≤ yrot ≤ tanα*xrot  &&  0 ≤ yrot ≤ cosα*b
end
function parallelogram_region(x, y, idx, dom::Domain)
    x0 = dom.region_params[idx][1]
    y0 = dom.region_params[idx][2]
    cosθ = dom.region_params[idx][3]
    sinθ = dom.region_params[idx][4]
    a = dom.region_params[idx][5]
    b = dom.region_params[idx][6]
    cosα = dom.region_params[idx][7]
    tanα = dom.region_params[idx][8]
    return parallelogram_region(x,y,x0,y0,cosθ,a,b,sinθ,cosα,tanα)
end
function parallelogram_region(x, y, idx, bnd::Boundary, sys::System)
    x0 = sys.domains[idx].domain_params[1]
    y0 = sys.domains[idx].domain_params[2]
    cosθ = sys.domains[idx].domain_params[3]
    sinθ = sys.domains[idx].domain_params[4]
    a = sys.domains[idx].domain_params[5]
    b = sys.domains[idx].domain_params[6]
    cosα = sys.domains[idx].domain_params[7]
    tanα = sys.domains[idx].domain_params[8]
    return parallelogram_region(x,y,x0,y0,cosθ,sinθ,a,b,cosα,tanα)
end


function rotate(x,y,cosθ,sinθ)
    if !iszero(sinθ)
        xrot = +cosθ*x - sinθ*y
        yrot = +sinθ*x + cosθ*y
    else
        xrot = x
        yrot = y
    end
    return xrot, yrot
end
