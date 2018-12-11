module Shapes

export Shape,
Circle,
Ellipse,
Square,
Rectangle,
Parallelogram


abstract type Shape end


"""
    Circle(R, x0, y0)
"""
struct Circle <: Shape
    R::Float64
    x0::Float64
    y0::Float64
    name::Symbol

    function Circle(R::Number,x0::Number,y0::Number)
        new(R,x0,y0, :circle)
    end
    function (c::Circle)(x,y)
        return hypot((x-c.x0), (y-c.y0)) ≤ c.R
    end
end


"""
    Ellipse(a, b, x0, y0, θ)
"""
struct Ellipse <: Shape
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64
    name::Symbol

    function Ellipse(a::Number,b::Number,x0::Number,y0::Number,θ::Number)
        new(float(a),float(b),float(x0),float(y0),float(θ),cos(θ),sin(θ), :ellipse)
    end
    function (e::Ellipse)(x,y)
        xrot, yrot = rotate(x-e.x0, y-e.y0, e.cosθ, e.sinθ)
        return hypot(xrot/e.a, yrot/e.b) ≤ 1
    end
end


"""
    Square(a,x0,y0,θ)
"""
struct Square <: Shape
    a::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64
    name::Symbol

    function Square(a::Number,x0::Number,y0::Number,θ::Number)
        new(float(a),float(x0),float(y0),float(θ),cos(θ),sin(θ), :square)
    end
    function (s::Square)(x,y)
        xrot, yrot = rotate(x-s.x0, y-s.y0, s.cosθ, s.sinθ)
        return 0 ≤ xrot ≤ a  &&  0 ≤ yrot ≤ a
    end
end


"""
    Rectangle(a, b, x0, y0, θ)
"""
struct Rectangle <: Shape
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64
    name::Symbol

    function Rectangle(a::Number,b::Number,x0::Number,y0::Number,θ::Number)
        new(float(a),float(b),float(x0),float(y0),float(θ),cos(θ),sin(θ), :rectangle)
    end

    function (r::Rectangle)(x,y)
        xrot, yrot = rotate(x-r.x0, y-r.y0, r.cosθ, r.sinθ)
        return 0 ≤ xrot ≤ r.a  &&  0 ≤ yrot ≤ r.a
    end
end


"""
    Parallelogram(a, b, α, x0, y0, θ)
"""
struct Parallelogram <: Shape
    a::Float64
    b::Float64
    α::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64
    tanα::Float64
    cosα::Float64
    name::Symbol

    function Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number)
        new(float(a),float(b),float(α),float(x0),float(y0),float(θ),cos(θ),sin(θ),tan(α),cos(α), :parallelogram)
    end

    function (p::Parallelogram)(x,y)
        xrot, yrot = rotate(x-p.x0, y-p.y0, p.cosθ, p.sinθ)
        return p.tanα*(xrot-p.a) ≤ yrot ≤ p.tanα*xrot  &&  0 ≤ yrot ≤ p.cosα*p.b
    end
end


"""
    rotate(x,y,cosθ,sinθ)
"""
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

end #module
