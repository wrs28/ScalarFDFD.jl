#TODO: fix simpson's rule boundary correction, for now I've just copied trapezoidal

################################################################################
## QUADRATURE
################################################################################
"""
    ∫f(z) dz = quadrature(sim, f)

multi-dimensional integral of f, uses trapezoidal rule if number of quadrature
points is even and simpson if odd in each dimension.
"""
function quadrature(sim::Simulation, f::AbstractArray{T,1}) where T
    if length(f)==prod(sim.dis.N_tr)
        quadrature(f, sim.dis.N_tr, sim.dis.dx, sim.bnd.bc)
    else
        quadrature(f, sim.dis.N, sim.dis.dx, sim.bnc.bc)
    end
end


"""
    ∫f(z) dz = quadrature(f::Array{T,1}, N, dx)

reshapes `f` before integrating.
"""
function quadrature(f::AbstractArray{T,1}, N::Array{Int,1}, dx, bnds::Array{Symbol,2}=fill(:none,2,2)) where T
    if ! (1 ∈ N)
        return quadrature(reshape(f, N[1], N[2]), dx, bnds)
    elseif N[1] == 1
        return quadrature(f[:], dx[2], bnds)
    else
        return quadrature(f[:], dx[1], bnds)
    end
end


"""
    ∫f(z) dz = quadrature(f::Array{T,N}, dx)
"""
function quadrature(f::AbstractArray{T,N}, dx::AbstractArray{Float64,1}, bnds::Array{Symbol,2}=fill(:none,2,2)) where T where N
    for i ∈ 1:N
        if isodd(size(f)[i])
            f = simpson(f,dx,i,bnds)
        else
            f = trapz(f,dx,i,bnds)
        end
    end
    return f[1]
end


function quadrature(f::AbstractArray{T,N}, dx::Float64, bnds) where T where N
    return quadrature(f, fill(dx, N), bnds)
end


function quadrature(f::AbstractArray{T,N}, dx::AbstractArray{Float64,1}, bnds::Symbol) where T where N
    return quadrature(f, dx, fill(bnds,2,2))
end


################################################################################
## INTEGRATION METHODS
################################################################################
"""
    ∫f(z) dz = trapz(f, dz, dim, bnds=[:none,:none])
"""
function trapz(f::AbstractArray{T,N}, dz::Array{Float64,1}, dim::Int, bnds::Array{Symbol,2}=fill(:none,2,2)) where T where N

    full_dims = [size(f)...]
    collapsed_dims = copy(full_dims)
    collapsed_dims[dim] = 1

    Q = zeros(T,collapsed_dims...)::Array{T,N}

    for i ∈ 2:full_dims[dim]-1
        Q[:] = Q[:] + 2selectdim(f,dim,i)[:]
    end
    Q[:] = Q[:] + selectdim(f,dim,1)[:]
    Q[:] = Q[:] + selectdim(f,dim,full_dims[dim])[:]

    for i ∈ 1:2
        i==1 ? idx=1 : idx=full_dims[dim]
        if bnds[i,dim] == :none
            nothing
        elseif bnds[i,dim] == :d
            Q[:] = Q[:] + selectdim(f,dim,idx)[:]/2
        elseif bnds[i,dim] ∈ [:n,:p]
            Q[:] = Q[:] + selectdim(f,dim,idx)[:]
        else
            throw(ArgumentError("trapz not defined for boundary condition $(sim.bnd.bc[i,dim])"))
        end
    end

    return Q*dz[dim]/2
end


"""
    ∫f(z) dz = simpson(f, dz, dim, bnds=[:none,:none])
"""
function simpson(f::AbstractArray{T,N}, dz::Array{Float64,1}, dim::Int, bnds::Array{Symbol,2}=fill(:none,2,2)) where T where N
    if iseven(size(f,dim))
        @warn "n=$(size(f,dim)) even, results not accurate."
    end

    full_dims = [size(f)...]
    collapsed_dims = copy(full_dims)
    collapsed_dims[dim] = 1

    Q = zeros(T,collapsed_dims...)::Array{T,N}

    for i ∈ 2:2:full_dims[dim]-1
        Q[:] = Q[:] + 4selectdim(f,dim,i)[:]
    end
    for i ∈ 3:2:full_dims[dim]-1
        Q[:] = Q[:] + 2selectdim(f,dim,i)[:]
    end
    Q[:] = Q[:] + selectdim(f,dim,1)[:]
    Q[:] = Q[:] + selectdim(f,dim,full_dims[dim])[:]

    for i ∈ 1:2
        i==1 ? idx=1 : idx=full_dims[dim]
        if bnds[i,dim] == :none
            nothing
        elseif bnds[i,dim] == :d
            Q[:] = Q[:] + 3selectdim(f,dim,idx)[:]/4
        elseif bnds[i,dim] ∈ [:n,:p]
            Q[:] = Q[:] + 3selectdim(f,dim,idx)[:]/2
        else
            throw(ArgumentError("simpson not defined for boundary condition $(sim.bnd.bc[i,dim])"))
        end
    end

    return Q*dz[dim]/3
end
