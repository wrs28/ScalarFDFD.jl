# TODO: check effect of polar coordinates on ihomogeneous systems.
# TODO: loop over size of M to apply BC's deeper in the interior in order to do photonic crystal waveguides as open boundaries
# TODO: make polar coordinates available in d=1, which is useful for, say, central-potential scattering

"""
    module DifferentialOperators

module for building gradient and laplacian in 1 or 2 dimensions subject to either periodic or robin boundary conditions.

For d=2, can specify either cartesian or polar coordinates.

Use `ezbc` to construct boundary conditions easily, or `RobinBoundaryCondition` and `PeriodicBoundaryCondition` for complete control.

`gradient` and `laplacian` give sparse matrices that implement the ∇ and ∇⋅(h∇) operators.
"""
module DifferentialOperators

using LinearAlgebra,
SparseArrays,
Bravais

export BoundaryCondition,
RobinBoundaryCondition,
PeriodicBoundaryCondition,
ezbc,
gradient,
laplacian,
isPolar,
isCartesian

# must define here for type stability later
function ⊗(A,B)
    kron(A,B)
end

################################################################################
### BOUNDARY CONDITIONS
################################################################################
"""
    BoundaryCondition

Abstract container for PeriodicBoundaryCondition and RobinBoundaryCondition

See also: [`RobinBoundaryCondition`](@ref), [`PeriodicBoundaryCondition`](@ref)
"""
abstract type BoundaryCondition end


"""
    pbc = PeriodicBoundaryCondition(N::Array{Int}, dimension, lattice=BravaisLattice())

periodic boundary conditions along `dimension` where the domain is periodic with `lattice`

Fields of `pbc` are:

* `lattice`, the Bravais lattice

* `N`, an array of the number of sites in each dimension

* `row_inds`, `col_inds`, `weights` give the location and coefficient of the connection between boundaries

* `shifts_a`, `shifts_b` give the number of dislocations along each of the Bravais vectors necessary to come back to the rectangular unit cell

* `dim` is the dimension along which the periodic condition applies

See also: [`RobinBoundaryCondition`](@ref), [`BoundaryCondition`](@ref),  [`ezbc`](@ref)
"""
struct PeriodicBoundaryCondition <: BoundaryCondition
    lattice::BravaisLattice
    N::Array{Int,1}
    row_inds::Array{Int,1}
    col_inds::Array{Int,1}
    weights::Array{Float64,1}
    shifts_a::Array{Int,1}
    shifts_b::Array{Int,1}
    dim::Int

    function PeriodicBoundaryCondition(N::Array{Int,1}, dim::Int, lattice::BravaisLattice=BravaisLattice(a = N[1], b = N[2]))

        if 0 ∈ N
            throw(ArgumentError("improperly initialized N = $N"))
        end

        a, α, v1 = lattice.a, lattice.α, lattice.v1
        b, β, v2 = lattice.b, lattice.β, lattice.v2

        width  = a*sin(β-α)/sin(β)
        height = b*sin(β)

        dx, dy = width/N[1], height/N[2]

        if (dim == 1 && isinf(a)) || (dim == 2 && isinf(b))
            return new(lattice, N, Array{Int}(undef,0), Array{Int}(undef,0), Array{Float64}(undef,0), Array{Int}(undef,0), Array{Int}(undef,0), dim)
        elseif !(dim ∈ [1,2])
            throw(ArgumentError("invalid dimensions $(dim)"))
        else
            nothing
        end

        if !isinf(a)
            start, stop = lattice.x0 + min(0, a*sin(β-α)/sin(β)), lattice.x0 + max(0, a*sin(β-α)/sin(β))
            x = LinRange(start + dx/2, stop - dx/2, N[1])
        else
            x = [lattice.x0]
        end

        if !isinf(b)
            start, stop = lattice.y0 + min(0, b*sin(β)), lattice.y0 + max(0, b*sin(β))
            y = LinRange(start + dy/2, stop - dy/2, N[2])
        else
            y = [lattice.y0]
        end
        if dim == 1
            p1, p2 = bravais_coordinates(x[1]-dx, y, lattice)
        else
            p1, p2 = bravais_coordinates(x, y[1]-dy, lattice)
        end

        Ma, Mb = -floor.(Int, p1/a), -floor.(Int, p2/b)
        if isinf(a)
            X = v1[1]*p1 + v2[1]*(p2 + Mb*b)
            Y = v1[2]*p1 + v2[2]*(p2 + Mb*b)
        elseif isinf(b)
            X = v1[1]*(p1 + Ma*a) + v2[1]*p2
            Y = v1[2]*(p1 + Ma*a) + v2[2]*p2
        else
            X = v1[1]*(p1 + Ma*a) + v2[1]*(p2 + Mb*b)
            Y = v1[2]*(p1 + Ma*a) + v2[2]*(p2 + Mb*b)
        end

        Ma += -floor.(Int, X/width)
        Mb += -floor.(Int, Y/height)

        x_inds1 = floor.(Int, X/dx .+ 1/2)
        x_inds2 = x_inds1 .+ 1

        y_inds1 = floor.(Int, Y/dy .+ 1/2)
        y_inds2 = y_inds1 .+ 1

        Cx1, Cx2 = abs.(X/dx .+ 1/2 - x_inds2), abs.(X/dx .+ 1/2 - x_inds1)
        cx1, cx2 = Cx1./(Cx1+Cx2), Cx2./(Cx1+Cx2)

        Cy1, Cy2 = abs.(Y/dy .+ 1/2 - y_inds2), abs.(Y/dy .+ 1/2 - y_inds1)
        cy1, cy2 = Cy1./(Cy1+Cy2), Cy2./(Cy1+Cy2)

        q, r, s = Array{Int}(undef,2), Array{Int}(undef,2), Array{Float64}(undef,2)
        t, u, v = Array{Int}(undef,2), Array{Int}(undef,2), Array{Float64}(undef,2)
        j, k, l = Array{Int}(undef,4), Array{Int}(undef,4), Array{Float64}(undef,4)

        I = Array{Int}(undef, 4N[mod1(dim+1,2)])
        J = Array{Int}(undef, 4N[mod1(dim+1,2)])
        V = Array{Float64}(undef, 4N[mod1(dim+1,2)])
        Na = Array{Int}(undef, 4N[mod1(dim+1,2)])
        Nb = Array{Int}(undef, 4N[mod1(dim+1,2)])

        for i ∈ 1:N[mod1(dim+1,2)]

            if dim == 1
                ind_x, ind_y = 1, i
            else
                ind_x, ind_y = i, 1
            end

            q[1], q[2] = ind_x, ind_x
            r[1], r[2] = mod1(x_inds1[i],N[1]), mod1(x_inds2[i],N[1])
            s[1], s[2] = cx1[i], cx2[i]

            t[1], t[2] = ind_y, ind_y
            u[1], u[2] = mod1(y_inds1[i],N[2]), mod1(y_inds2[i],N[2])
            v[1], v[2] = cy1[i], cy2[i]

            if 1 ∈ N
                j[1:2], k[1:2], l[1:2] = findnz( sparse(t, u, v, N[2], N[2]) ⊗ sparse(q, r, s, N[1], N[1]) )
                j[3:4] .= k[3:4] .= 1
                l[3:4] .= 0
            else
                j[:], k[:], l[:] = findnz( sparse(t, u, v, N[2], N[2]) ⊗ sparse(q, r, s, N[1], N[1]) )
            end

            I[(4(i-1)+1):(4(i-1)+4)] = j
            J[(4(i-1)+1):(4(i-1)+4)] = k
            V[(4(i-1)+1):(4(i-1)+4)] = l
            Na[(4(i-1)+1):(4(i-1)+4)] .= Ma[i]
            Nb[(4(i-1)+1):(4(i-1)+4)] .= Mb[i]
        end

        return new(lattice, N, I, J, V, Na, Nb, dim)
    end
    function PeriodicBoundaryCondition(N::Int, dim::Int, lattice::BravaisLattice=BravaisLattice(a = dim==1 ? float(N) : Inf, b = dim==2 ? float(N) : Inf))
        if dim == 1
            N = [N,1]
        else
            N = [1,N]
        end
        return PeriodicBoundaryCondition(N, dim, lattice)
    end
end # struct PeriodicBoundaryCondition


"""
    rbc = RobinBoundaryCondition(;α,β,p=[0.,0],q=[0.,0],g=[0.,0])

Robin Boundary conditions A*ϕ+B*∇ϕ = G, defaulting to Dirichlet in 1 dim

A = [α[1] p[1]      B = [β[1] q[1]      G = [g[1]
     p[2] α[2]]          q[2] β[2]]          g[2]]
"""
struct RobinBoundaryCondition{T} <: BoundaryCondition
    α::Array{Array{T,2},1}
    β::Array{Array{T,2},1}
    p::Array{Array{T,2},1}
    q::Array{Array{T,2},1}
    g::Array{Array{T,1},1}
    A::Array{T,2}
    B::Array{T,2}
    G::Array{T,1}

    function RobinBoundaryCondition(
        α::Array{Array{T,2},1},
        β::Array{Array{T,2},1},
        p::Array{Array{T,2},1},
        q::Array{Array{T,2},1},
        g::Array{Array{T,1},1}
        ) where T <: Number

        A = [α[1] p[1]
             p[2] α[2]]

        B = [β[1] q[1]
             q[2] β[2]]

        G = vcat(g[1], g[2])

        return new{T}(α, β, p, q, g, A, B, G)
    end
    function RobinBoundaryCondition(α::Array{T,1}, β::Array{Array{T,2},1}, p::Array{Array{T,2},1}, q::Array{Array{T,2},1}, g::Array{Array{T,1},1}) where T<:Number
        return RobinBoundaryCondition([hcat(α[i],) for i ∈ eachindex(α)], β, p, q, g)
    end
    function RobinBoundaryCondition(α, β::Array{T,1}, p::Array{Array{T,2},1}, q::Array{Array{T,2},1}, g::Array{Array{T,1},1}) where T<:Number
        return RobinBoundaryCondition(α, [hcat(β[i],) for i ∈ eachindex(β)], p, q, g)
    end
    function RobinBoundaryCondition(α, β, p::Array{T,1}, q::Array{Array{T,2},1}, g::Array{Array{T,1},1}) where T<:Number
        return RobinBoundaryCondition(α, β, [hcat(p[i],) for i ∈ eachindex(p)], q, g)
    end
    function RobinBoundaryCondition(α, β, p, q::Array{T,1}, g::Array{Array{T,1},1}) where T<:Number
        return RobinBoundaryCondition(α, β, p, [hcat(q[i],) for i ∈ eachindex(q)], g)
    end
    function RobinBoundaryCondition(α, β, p, q, g::Array{T,1}) where T<:Number
        return RobinBoundaryCondition(α, β, p, q, [[g[i]] for i ∈ eachindex(g)])
    end
    function RobinBoundaryCondition(;α,β,
        p=zeros(Float64,size(α)),
        q=zeros(Float64,size(α)),
        g=zeros(Float64,size(α)))
        return RobinBoundaryCondition(α,β,p,q,g)
    end
    function RobinBoundaryCondition(bc::Array{Symbol,1})
        α = Array{Float64}(undef,2)
        β = Array{Float64}(undef,2)
        p = zeros(Float64,2)
        q = zeros(Float64,2)
        g = zeros(Float64,2)
        for j ∈ 1:2
            if bc[j] == :n
                α[j] = 0.
                β[j] = 1.
            elseif bc[j] == :d
                α[j] = 1.
                β[j] = 0.
            else
                throw(ArgumentError("unrecognized boundary condition $bc"))
            end
        end
        return RobinBoundaryCondition(α, β, p, q, g)
    end
    function RobinBoundaryCondition(bc1::Symbol,bc2::Symbol)
        return RobinBoundaryCondition([bc1,bc2])
    end
    function RobinBoundaryCondition(bc::Symbol)
        return RobinBoundaryCondition([bc,bc])
    end
end # struct RobinBoundaryCondition


"""
    bcs = ezbc(bc1, bc2, N=[0,0])
    bcs = ezbc(bcs, N=[0,0])

wrapper for easy construction of boundary conditions, `bc1` for dimension 1, `bc2` for dimension 2.

`bc1`, `bc2` can by specified as symbols or arrays of symbols. If symbol, it applies to both boundaries in that dimension.

See also: [`RobinBoundaryCondition`](@ref), [`PeriodicBoundaryCondition`](@ref)
"""
function ezbc(bc1::Array{Symbol,1},bc2::Array{Symbol,1}, N::Array{Int}=[1,1])
    if issubset(bc1,[:n, :d])
        BC1 = RobinBoundaryCondition(bc1[1],bc1[2])
    elseif bc1 == [:p,:p]
        BC1 = PeriodicBoundaryCondition(N, 1)
    else
        throw(ArgumentError("unrecognized boundary specification $bc1"))
    end
    if issubset(bc2,[:n, :d])
        BC2 = RobinBoundaryCondition(bc2[1],bc2[2])
    elseif bc2 == [:p,:p]
        BC2 = PeriodicBoundaryCondition(N, 2)
    else
        throw(ArgumentError("unrecognized boundary specification $bc2"))
    end
    return (BC1,BC2)
end
function ezbc(bc1::Symbol,bc2::Array{Symbol,1}, N::Array{Int}=[0,0])
    return ezbc([bc1, bc1],bc2, N)
end
function ezbc(bc1::Array{Symbol,1},bc2::Symbol, N::Array{Int}=[0,0])
    return ezbc(bc1, [bc2,bc2], N)
end
function ezbc(bc1::Symbol,bc2::Symbol,N::Array{Int})
    return ezbc([bc1, bc1], [bc2, bc2], N)
end

function ezbc(bc::Array{Symbol,1}, N::Int=0)
    if issubset(bc,[:n, :d])
        return RobinBoundaryCondition(bc)
    elseif issubset(bc, [:p, :p])
        return PeriodicBoundaryCondition(N, findfirst(N.!==1))
    else
        throw(ArgumentError("unrecognized boundary specification $bc"))
    end
end
function ezbc(bc1::Symbol,bc2::Symbol, N::Int)
    return ezbc([bc1,bc2], N)
end
function ezbc(bc::Symbol, N::Int=0)
    return ezbc([bc,bc], N)
end


################################################################################
####### GRADIENTS
################################################################################
"""
    ∇ = grad(N::Int, dx, bc; polarity=:central, ka=0, kb=0, coordinate_system=:cart)

`bc` can be an instance of `RobinBoundaryCondition` or `PeriodicBoundaryCondition`.

`polarity` is one of `:forward`, `:backward`, or `:central`

`N` is number of lattice sites, `dx` is lattice spacing

See also: [`ezbc`](@ref), [`laplacian`](@ref)
"""
function grad(N::Int, dx::Real, bc::T; polarity::Symbol=:central, ka::Number=0, kb::Number=0, coordinate_system::Symbol=:cart) where T<:BoundaryCondition
    ∇ = grad_sans_bc(N, dx, polarity, coordinate_system)
    ∇, S = grad_apply_bc(∇, [N,1], dx[1], bc, 1, polarity, ka, kb, coordinate_system)
    return ∇, S
end # 1d grad


"""
    ∇₁, ∇₂ =  grad(N::Array{Int}, dx::Array{Float64}, bcs; polarity=:central, ka=0, kb=0, coordinate_system=:cart)

2-dim gradients with `N[1]`, `N[2]` points, lattice spacing `dx[1], dx[2]`

See also: [`ezbc`](@ref), [`laplacian`](@ref)
"""
function grad(N::Array{Int}, dx::Array{T}, bcs::Tuple{U,V};
                polarity::Symbol=:central, ka::Number=0, kb::Number=0, coordinate_system::Symbol=:cart) where T<:Real where U<:BoundaryCondition where V<:BoundaryCondition

    ∇1 = grad_sans_bc(N[1], dx[1], polarity, coordinate_system)
    ∇2 = grad_sans_bc(N[2], dx[2], polarity, coordinate_system)

    𝕀1, 𝕀2 = sparse(I, N[1], N[1]), sparse(I, N[2], N[2])
    ∇₁ = 𝕀2 ⊗ ∇1
    if isPolar(coordinate_system)
        ∇₂ = ∇2⊗spdiagm(0 => 1 ./(dx[2]*(1/2 .+ (0:N[1]-1))))
    elseif isCartesian(coordinate_system)
        ∇₂ = ∇2 ⊗ 𝕀1
    else
        throw(ArgumentError("unrecognized coordinate system $coordinate_system"))
    end

    ∇₁, S₁ = grad_apply_bc(∇₁, N, dx, bcs, 1, polarity, ka, kb, coordinate_system)
    ∇₂, S₂ = grad_apply_bc(∇₂, N, dx, bcs, 2, polarity, ka, kb, coordinate_system)

    return ∇₁, ∇₂, S₁, S₂
end # 2d grad


"""
    grad_sans_bc(N, dx, polarity, coordinate_system)
"""
function grad_sans_bc(N::Int, dx::Real, polarity::Symbol, coordinate_system::Symbol)
    if isCentral(polarity)
        dx = 2dx
        I₁, J₁ = Array(2:N), Array(1:N-1)
        V₁ = fill(-1/dx, N-1)
        I₂, J₂ = Array(1:N-1), Array(2:N)
        V₂ = fill(+1/dx, N-1)
    elseif isForward(polarity)
        I₁, J₁ = Array(1:N), Array(1:N)
        V₁ = fill(-1/dx, N)
        I₂, J₂ = Array(1:N-1), Array(2:N)
        V₂ = fill(+1/dx, N-1)
    elseif isBackward(polarity)
        I₁, J₁ = Array(2:N), Array(1:N-1)
        V₁ = fill(-1/dx, N-1)
        I₂, J₂ = Array(1:N), Array(1:N)
        V₂ = fill(+1/dx, N)
    end
    ∇ = sparse(vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N, N, +)
    return ∇
end # grad_sans_bc


"""
    grad_apply_bc(∇, N, dx, bc, dim, polarity, ka, kb, coordinate_system)
"""
function grad_apply_bc(∇, N, dx, bc, dim::Int, polarity::Symbol, ka::Number, kb::Number, coordinate_system::Symbol)
    return apply_bc(∇, :gradient, N, dx, bc, dim, 1, polarity, ka, kb, coordinate_system)
end # grad_apply_bc


################################################################################
####### LAPLACIANS
################################################################################
"""
    ∇², S = laplacian(N::Array{Int}, dx, bcs::Tuple; ka=0, kb=0, coordinate_system=:cart, h=1)

Compute ∇⋅(`h`∇) on a 1-dim lattice with `N` sites, spacing `dx`, subject to boundary conditions contained in `bc`.

`ka`, `kb` are Floquet phases.

`S` is still experimental, should specify inhomogneous boundary term if there is one specified in `bcs`

See also: [`ezbc`](@ref), [`grad`](@ref)
"""
function laplacian(N::Int, dx, bc::BoundaryCondition; ka::Number=0, kb::Number=0, h=1)
    ∇² = laplacian_sans_bc(N[1], dx[1], h, :cart)
    ∇², S = laplacian_apply_bc(∇², [N[1],1], dx[1], bc, 1, h, ka, kb, :cart)
    return ∇², S
end # 1d laplacian


"""
    ∇², S = laplacian(N::Array{Int}, dx, bcs::Tuple; ka=0, kb=0, coordinate_system=:cart, h=1)

Compute ∇⋅(`h`∇) on a 2-dim lattice with `N[1]`×`N[2]` sites, spacings `dx[1]`, `dx[2]`, subject to boundary conditions contained in `bcs`.

`ka`, `kb` are Floquet phases, `coordinate_system` specifies cartesian or polar.

`S` is still experimental, should specify inhomogneous boundary term if there is one specified in `bcs`

See also: [`ezbc`](@ref), [`grad`](@ref)
"""
function laplacian(N::Array{Int}, dx, bcs::Tuple{U, V};
            ka::Number=0, kb::Number=0, coordinate_system::Symbol=:cart, h=1) where U<:BoundaryCondition where V<:BoundaryCondition

    ∇1² = laplacian_sans_bc(N[1], dx[1], h, coordinate_system)
    ∇2² = laplacian_sans_bc(N[2], dx[2], h, :cart)
    𝕀1, 𝕀2 = sparse(I, N[1], N[1]), sparse(I, N[2], N[2])
    ∇₁², ∇₂² = 𝕀2 ⊗ ∇1², ∇2² ⊗ 𝕀1

    ∇₁², S1 = laplacian_apply_bc(∇₁², N, dx, bcs, 1, h, ka, kb, coordinate_system)
    ∇₂², S2 = laplacian_apply_bc(∇₂², N, dx, bcs, 2, h, ka, kb, coordinate_system)

    if isPolar(coordinate_system)
        r⁻² = 𝕀2⊗sparse(1:N[1],1:N[1], 1 ./(dx[1]*(1/2 .+ (0:N[1]-1))).^2, N[1], N[1])
        ∇₂² = r⁻²*∇₂²
    end

    return ∇₁² + ∇₂², S1+S2
end
function laplacian(N::Array{Int}, dx::Real, bcs::Tuple{U, V};
            ka::Number=0, kb::Number=0, coordinate_system::Symbol=:cart, h=1) where U<:BoundaryCondition where V<:BoundaryCondition
    ∇², S = laplacian(N, [dx, dx], bcs; coordinate_system=coordinate_system, h=h, ka=ka, kb=kb)
    return ∇², S
end # 2d laplacian


"""
    ∇h∇ = laplacian_sans_bc(N, dx, h, coordinate_system)
"""
function laplacian_sans_bc(N::Int, dx::Real, h::AbstractArray, coordinate_system::Symbol)
    I₁, J₁ = Array(2:N), Array(1:N-1)
    V₁ = +1h[2:end-1]/dx^2
    if isPolar(coordinate_system)
        V₁ = ((I₁.-1)./sqrt.((I₁.-1).^2 .- 1/4)).*V₁
    end

    I₂, J₂ = Array(1:N), Array(1:N)
    V₂ = -(h[1:end-1]+h[2:end])/dx^2

    I₃, J₃, V₃ = J₁, I₁, V₁

    ∇² = sparse(vcat(I₁,I₂,I₃), vcat(J₁,J₂,J₃), vcat(V₁,V₂,V₃), N, N, +)
    return ∇²
end
function laplacian_sans_bc(N::Int, dx::Real, h::Number, coordinate_system::Symbol)
    return laplacian_sans_bc(N, dx, fill(h,N+1), coordinate_system)
end # laplacian_sans_bc


"""
    ∇², S = laplacian_apply_bc(∇², N, dx, bc, dim, h, ka, kb, coordinate_system)
"""
function laplacian_apply_bc(∇², N, dx, bc::Tuple{T,U}, dim, h, ka, kb, coordinate_system) where T<:BoundaryCondition where U<:BoundaryCondition
    ∇², S = apply_bc(∇², :laplacian, N, dx, bc, dim, h, :central, ka, kb, coordinate_system)
    return ∇², S
end # laplacian_apply_bc


########################################################################################
### APPLICATION OF BOUNDARY CONDITIONS
########################################################################################
"""
    apply_bc(D, operator::Symbol, N::Array, dx::Real, bc<:BoundaryCondition, dim, h, polarity, ka, kb, coordinate_system)
"""
function apply_bc(D, operator::Symbol, N::Array{Int,1}, dx::Real, bc::RobinBoundaryCondition, dim::Int, h::AbstractArray{V,1},
            polarity::Symbol, ka::Number, kb::Number, coordinate_system::Symbol) where V<:Number

    if !all(size(bc.α[1]) .== fill(N[mod1(dim+1,2)],2))
        A = diagm(0=>vcat(fill(bc.α[1][1],N[mod1(dim+1,2)]),fill(bc.α[2][1],N[mod1(dim+1,2)])))
    else
        A = bc.A
    end
    if !all(size(bc.β[1]) .== fill(N[mod1(dim+1,2)],2) )
        B = diagm(0=>vcat(fill(bc.β[1][1],N[mod1(dim+1,2)]),fill(bc.β[2][1],N[mod1(dim+1,2)])))
    else
        B = bc.B
    end
    if length(bc.g[1]) !== N[mod1(dim+1,2)]
        G = vcat(fill(bc.g[1][1],N[mod1(dim+1,2)]),fill(bc.g[2][1],N[mod1(dim+1,2)]))
    else
        G = bc.G
    end
    C = inv(A+B/dx)
    M = -C*(A-B/dx)
    S = +2C*G

    if operator==:laplacian
        dx = dx^2
        sgn = +1
        if isPolar(coordinate_system)
            multiple_inner = 0
            multiple_outer = (N[1]-1)/sqrt((N[1]-1)^2-1/4)
            multiple_outer_g = sqrt(dx*(N[1]+1/2))
        else
            multiple_inner=1
            multiple_outer=1
            multiple_outer_g=1
        end
    elseif operator==:gradient
        dx = isCentral(polarity) ? 2dx : dx
        sgn = -1
        multiple_inner=1
        multiple_outer=1
        multiple_outer_g=1
    else
        throw(ArgumentError("unrecognized differential operator $operator"))
    end

    if dim==1
        BC  = multiple_inner*sgn*h[1]*M[1:end÷2,1:end÷2]⊗sparse([1], [1], [1/dx], N[1], N[1])
        BC += sgn*h[1]*M[end÷2+1:end,1:end÷2]⊗sparse([1], [N[1]], [1/dx], N[1], N[1])
        BC += +h[end]*M[1:end÷2,end÷2+1:end]⊗sparse([N[1]], [1], [1/dx], N[1], N[1])
        BC += +multiple_outer*h[end]*M[end÷2+1:end,end÷2+1:end]⊗sparse([N[1]], [N[1]], [1/dx], N[1], N[1])
        BG  = multiple_inner*S[1:end÷2]⊗vcat(sgn/dx,zeros(N[1]-1))
        BG += multiple_outer*multiple_outer_g*S[end÷2+1:end]⊗vcat(zeros(N[1]-1),1/dx)
    elseif dim==2
        BC  = sgn*h[1]*sparse([1], [1], [1/dx], N[2], N[2])⊗M[1:end÷2,1:end÷2]
        BC += sgn*h[1]*sparse([1], [N[2]], [1/dx], N[2], N[2])⊗M[end÷2+1:end,1:end÷2]
        BC += +h[end]*sparse([N[2]], [1], [1/dx], N[2], N[2])⊗M[1:end÷2,end÷2+1:end]
        BC += +h[end]*sparse([N[2]], [N[2]], [1/dx], N[2], N[2])⊗M[end÷2+1:end,end÷2+1:end]
        BG  = vcat(sgn/dx,zeros(N[2]-1))⊗S[1:end÷2]
        BG += vcat(zeros(N[2]-1),1/dx)⊗S[end÷2+1:end]
    else
        throw(ArgumentError("unrecognized dimension $dim"))
    end
    return D + BC, -BG
end # apply_bc RobinBoundaryCondition

function apply_bc(D, operator::Symbol, N::Array{Int}, dx::Real, bc::PeriodicBoundaryCondition, dim::Int, h::AbstractArray{T,1},
            polarity::Symbol, ka::Number, kb::Number, coordinate_system::Symbol) where T<:Number

    if operator==:laplacian
        dx = dx^2
        sgn = +1
    elseif operator==:gradient
        dx = isCentral(polarity) ? 2dx : dx
        sgn = -1
    else
        throw(ArgumentError("unrecognized differential operator $operator"))
    end

    N, lattice, shifts_a, shifts_b = bc.N, bc.lattice, bc.shifts_a, bc.shifts_b
    rows, cols, weights = bc.row_inds, bc.col_inds, bc.weights

    if !isinf(lattice.a) && !isinf(lattice.b)
        ϕ = -shifts_a*ka*lattice.a - shifts_b*kb*lattice.b
    elseif !isinf(lattice.b)
        ϕ = -shifts_b*kb*lattice.b
    elseif !isinf(lattice.a)
        ϕ = -shifts_a*ka*lattice.a
    else
        ϕ = 0
    end
    BC = sparse(vcat(rows,cols), vcat(cols,rows), vcat(sgn*weights.*exp.(+1im*ϕ)/dx,+weights.*exp.(-1im*ϕ)/dx), prod(N), prod(N), +)
    return D + BC, zeros(Float64,size(D,1))
end # apply_bc PeriodicBoundaryCondition

"""
    apply_bc(D, operator::Symbol, N::Array, dx::Array, bc::Tuple, dim, h, polarity, ka, kb, coordinate_system)
"""
function apply_bc(D, operator, N, dx, bc, dim, h::Number, polarity::Symbol, ka::Number, kb::Number, coordinate_system::Symbol)
    return apply_bc(D, operator, N, dx, bc, dim, fill(h,N[dim]), polarity, ka, kb, coordinate_system)
end
function apply_bc(D, operator::Symbol, N::Array{Int}, dx::Array{U,1}, bc::Tuple{T,S}, dim::Int, h::AbstractArray{V,1},
            polarity::Symbol, ka::Number, kb::Number, coordinate_system::Symbol) where U<:Real where V<:Real where T<:BoundaryCondition where S<:BoundaryCondition
    return apply_bc(D, operator, N, dx[dim], bc[dim], dim, h[dim], polarity, ka, kb, coordinate_system)
end # apply_bc


################################################################################################
### ISS
################################################################################################
"""
    bool = isCartesian(coordinate_system::Symbol)
"""
function isCartesian(coordinate_system::Symbol)
    validNames = [:Cartesian, :cartesian, :Cart, :cart]
    return coordinate_system ∈ validNames
end


"""
    bool = isPolar(coordinate_system::Symbol)
"""
function isPolar(coordinate_system::Symbol)
    validNames = [:Polar, :polar, :Pol, :pol, :Cylindrical, :cylindrical, :Cyl, :cyl]
    return coordinate_system ∈ validNames
end


"""
    bool = isCentral(polarity::Symbol)
"""
function isCentral(polarity::Symbol)
    validNames = [:central, :Central, :center, :Center, :symmetric, :Symmetric, :cen, :Cen, :sym, :Sym, :c, :s, :C, :S]
    return polarity ∈ validNames
end


"""
    bool = isBackward(polarity::Symbol)
"""
function isBackward(polarity::Symbol)
    validNames = [:backward, :Backward, :back, :Back, :b, :B, :backwards, :Backwards]
    return polarity ∈ validNames
end

"""
    bool = isForward(polarity::Symbol)
"""
function isForward(polarity::Symbol)
    validNames = [:forward, :Forward, :for, :For, :f, :F, :forwards, :Forwards]
    return polarity ∈ validNames
end

end # module
